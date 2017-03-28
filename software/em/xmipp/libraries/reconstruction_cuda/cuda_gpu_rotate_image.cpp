
//Host includes
#include "cuda_gpu_rotate_image.h"
#include <iostream>
#include <stdio.h>
//CUDA includes
#include <cuda_runtime.h>
#include "cuda_copy_data.cpp"


// 2D float texture
texture<float, cudaTextureType2D, cudaReadModeElementType> texRef;

// 3D float texture
texture<float, cudaTextureType3D, cudaReadModeElementType> texRefVol;


//CUDA functions

// Cubic B-spline function
// The 3rd order Maximal Order and Minimum Support function, that it is maximally differentiable.
__device__ float bspline(float t)
{
	t = fabs(t);
	const float a = 2.0f - t;

	if (t < 1.0f) return 2.0f/3.0f - 0.5f*t*t*a;
	else if (t < 2.0f) return a*a*a / 6.0f;
	else return 0.0f;
}


//! Bicubic interpolated texture lookup, using unnormalized coordinates.
//! Straight forward implementation, using 16 nearest neighbour lookups.
//! @param tex  2D texture
//! @param x  unnormalized x texture coordinate
//! @param y  unnormalized y texture coordinate
__device__ float cubicTex2DSimple(texture<float, cudaTextureType2D, cudaReadModeElementType> tex, float x, float y)
{
	// transform the coordinate from [0,extent] to [-0.5, extent-0.5]
	const float2 coord_grid = make_float2(x - 0.5f, y - 0.5f);
	float2 index = make_float2(floor(coord_grid.x), floor(coord_grid.y));
	const float2 fraction = make_float2(coord_grid.x - index.x, coord_grid.y - index.y);
	index.x += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]
	index.y += 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]

	float result = 0.0f;
	for (float y=-1; y < 2.5f; y++)
	{
		float bsplineY = bspline(y-fraction.y);
		float v = index.y + y;
		for (float x=-1; x < 2.5f; x++)
		{
			float bsplineXY = bspline(x-fraction.x) * bsplineY;
			float u = index.x + x;
			result += bsplineXY * tex2D(tex, u, v);
		}
	}
	return result;
}


__global__ void
rotate_kernel_normalized(float *output, size_t Xdim, size_t Ydim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    u -= 0.5f;
    v -= 0.5f;

    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = tex2D(texRef, tu, tv);

}


__global__ void
rotate_kernel_unnormalized(float *output, size_t Xdim, size_t Ydim, float ang)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int y = blockDim.y * blockIdx.y + threadIdx.y;

    // Transform coordinates
    float u = x / (float)Xdim;
    float v = y / (float)Ydim;
    u -= 0.5f;
    v -= 0.5f;

    float tu = u * cosf(ang) - v * sinf(ang) + 0.5f;
    float tv = v * cosf(ang) + u * sinf(ang) + 0.5f;

    tu = tu*(float)Xdim;
    tv = tv*(float)Ydim;

    // Read from texture and write to global memory
   	output[y * Xdim + x] = cubicTex2DSimple(texRef, tu, tv);
}



void cuda_rotate_image(float *image, float *rotated_image, size_t Xdim, size_t Ydim, size_t Zdim, float ang, int interp){

	std::cerr  << "Inside CUDA function " << ang << std::endl;

	//CUDA code
	size_t matSize=Xdim*Ydim*Zdim*sizeof(float);
	struct cudaPitchedPtr bsplineCoeffs, cudaOutput;

	bsplineCoeffs = CopyVolumeHostToDevice(image, (uint)Xdim, (uint)Ydim, Zdim);

	if(Zdim==1){

		// Init texture
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
		cudaArray* cuArray;
		cudaMallocArray(&cuArray, &channelDesc, Xdim, Ydim);
		// Copy to device memory some data located at address h_data in host memory
		cudaMemcpy2DToArray(cuArray, 0, 0, bsplineCoeffs.ptr, bsplineCoeffs.pitch, Xdim * sizeof(float), Ydim, cudaMemcpyDeviceToDevice);

		// Bind the array to the texture reference
		cudaBindTextureToArray(texRef, cuArray, channelDesc);

	}/*else if (Zdim>1){

    	cudaMalloc3DArray(&cuArray, &channelDesc, make_cudaExtent(Xdim*sizeof(float),Ydim,Zdim), 0);
    	cudaMemcpy3DParms p = {0};
    	p.extent   = make_cudaExtent(Xdim*sizeof(float),Ydim,Zdim);
    	p.srcPtr   = make_cudaPitchedPtr((void*)image, Xdim * sizeof(float), Ydim, Zdim);
    	p.dstArray = *cuArray;
    	p.kind     = cudaMemcpyHostToDevice;
    	cudaMemcpy3D(&p);
    	// bind array to 3D texture
    	cudaBindTextureToArray(texRefVol, cuArray, channelDesc);

    }*/

    // Specify texture object parameters
    texRef.addressMode[0] = cudaAddressModeWrap;
    texRef.addressMode[1] = cudaAddressModeWrap;
    if(Zdim>1){
    	texRef.addressMode[2] = cudaAddressModeWrap;
    }
    if (interp==0){
    	texRef.filterMode = cudaFilterModePoint;
    }else{
    	texRef.filterMode = cudaFilterModeLinear;
    }
    if (interp<2){
    	texRef.normalized = true;
    }else{
    	texRef.normalized = false;
    }



    // Allocate result of transformation in device memory
    float *d_output;
	cudaMalloc((void **)&d_output, matSize);


	//Kernel
	int numTh = 32;
	const dim3 blockSize(numTh, numTh, 1);
	int numBlkx = (int)(Xdim)/numTh;
	if((Xdim)%numTh>0){
		numBlkx++;
	}
	int numBlky = (int)(Ydim)/numTh;
	if((Ydim)%numTh>0){
	  numBlky++;
	}
	const dim3 gridSize(numBlkx, numBlky, 1);

	if(interp<2){
		rotate_kernel_normalized<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, ang);
	}else{
		rotate_kernel_unnormalized<<<gridSize, blockSize>>>(d_output, Xdim, Ydim, ang);
	}

	cudaDeviceSynchronize();

	cudaMemcpy(rotated_image, d_output, matSize, cudaMemcpyDeviceToHost);


	cudaFree(cuArray);
	cudaFree(d_output);

}
