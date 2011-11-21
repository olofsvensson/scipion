package xmippij;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileInfo;
import ij.process.ImageProcessor;
import ij.process.ByteProcessor;
import ij.process.ShortProcessor;
import ij.process.FloatProcessor;
import ij.process.StackStatistics;
import java.io.File;
import java.util.LinkedList;
import xmipp.ImageGeneric_;
import xmipp.Projection;
import xmipp.MDLabel;
import xmipp.MetaData;


/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author Juanjo Vega
 */
public class XmippImageConverter_ {

    public static ImagePlus convertToImageJ(String filename) throws Exception {
        ImageGeneric_ image = new ImageGeneric_(filename);
        return convertToImageJ(image);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image) throws Exception {
        return convertToImageJ(image, ImageGeneric_.ALL_SLICES);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image, int nslice) throws Exception {
        return convertToImageJ(image, image.getXDim(), image.getYDim(), nslice);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image, long nimage) throws Exception {
        return convertToImageJ(image, image.getXDim(), image.getYDim(), nimage);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image, int width, int height) throws Exception {
        return convertToImageJ(image, width, height, ImageGeneric_.ALL_SLICES, ImageGeneric_.ALL_IMAGES);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image, int width, int height, int nslice) throws Exception {
        return convertToImageJ(image, width, height, nslice, ImageGeneric_.ALL_IMAGES);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image, int width, int height, long nimage) throws Exception {
        return convertToImageJ(image, width, height, ImageGeneric_.ALL_SLICES, nimage);
    }

    public static ImagePlus convertToImageJ(ImageGeneric_ image, int width, int height, int nslice, long nimage) throws Exception {
        ImageStack is = new ImageStack(width, height);
        ProcessorCreator_ pc = null;
        switch (image.getDataType()) {
            case ImageGeneric_.Float:
            case ImageGeneric_.Double:
                pc = new FloatProcessorCreator_();
                break;
            case ImageGeneric_.SChar:
            case ImageGeneric_.UChar:
                pc = new ByteProcessorCreator_();
                break;
            case ImageGeneric_.Short:
            case ImageGeneric_.UShort:
                pc = new ShortProcessorCreator_();
                break;
            default:
                pc = new FloatProcessorCreator_();
        }

        boolean retrieveAllImages = nimage == ImageGeneric_.ALL_IMAGES;
        boolean retrieveAllSlices = nslice == ImageGeneric_.ALL_SLICES;
        long n = retrieveAllImages ? ImageGeneric_.FIRST_IMAGE : nimage;

        for (; n <= image.getNDim(); n++) {
            System.out.println(" +++ Reading image: " + n);
            image.read(width, height, n);

            if (image.isPSD()) {
                image.convertPSD(image.getUseLogarithm());
            }

            // Read volume.
            int slice = retrieveAllSlices ? ImageGeneric_.FIRST_SLICE : nslice;
            for (; slice <= image.getZDim(); slice++) {
                System.out.println("\t /// Reading slice: " + slice);
                ImageProcessor processor = pc.getProcessor(image, slice);
                is.addSlice("", processor);

                // If just one image, breaks loop.
                if (!retrieveAllSlices) {
                    break;
                }
            }

            // If just one image, breaks loop.
            if (!retrieveAllImages) {
                break;
            }
        }

        return getNormalizedImagePlus(image.getFilename(), is);
    }

    public static void main(String args[]) {
        try {
            //String filename = "/home/jvega/Escritorio/imgs_Roberto/kk.stk"; // STACK
            //String filename = "/home/jvega/Escritorio/imgs_Roberto/kk.vol";   // VOLUME
            String filename = "/home/jvega/xmipp_ctf.psd"; // PSD
            ImageGeneric_ image = new ImageGeneric_(filename);

            convertToImageJ(image).show();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public static ImagePlus convertToImageJ(MetaData md) throws Exception {
        LinkedList<String> missing = new LinkedList<String>();
        ImagePlus imp = null;

        if (md.containsLabel(MDLabel.MDL_IMAGE)) {
            ImageStack is = null;

            long ids[] = md.findObjects();

            for (long id : ids) {
                String filename = md.getValueString(MDLabel.MDL_IMAGE, id, true);

                try {
                    ImagePlus slice = convertToImageJ(filename);

                    if (is == null) {
                        is = new ImageStack(slice.getWidth(), slice.getHeight());
                    }

                    is.addSlice(filename, slice.getProcessor());
                } catch (Exception ex) {
                    ex.printStackTrace(System.err);
                    missing.add(filename);
                }
            }

            imp = getNormalizedImagePlus(md.getFilename(), is);
        }

        // Tells user about missing files.
        if (!missing.isEmpty()) {
            String message = "There are missing files:\n";
            for (int i = 0; i < missing.size(); i++) {
                message += missing.get(i) + "\n";
            }

            IJ.error(message);
        }

        return imp;
    }

    public static ImagePlus convertToImageJ(Projection projection, String title) {
        FloatProcessor processor = new FloatProcessor(projection.getXsize(), projection.getYsize(), projection.getData());
        return new ImagePlus(title, processor);
    }

    public static ImagePlus getNormalizedImagePlus(String filename, ImageStack is) {
        ImagePlus imp = new ImagePlus(filename, is);
        // Sets associated file info.
        File f = new File(filename);
        FileInfo fi = new FileInfo();
        fi.directory = f.getParent();
        fi.fileName = f.getName();
        imp.setFileInfo(fi);

        // Normalize by default
        StackStatistics ss = new StackStatistics(imp);
        imp.getProcessor().setMinAndMax(ss.min, ss.max);
        return imp;
    }

    public static void revert(ImagePlus imp, String path) throws Exception {
        ImageGeneric_ image = new ImageGeneric_(path);
        ImagePlus imp2 = XmippImageConverter_.convertToImageJ(image);

        imp.setStack(imp.getTitle(), imp2.getImageStack());
//        int w = image.xSize;
//        int h = image.ySize;
//        int d = image.zSize;
//        long n = image.nSize;
//
//        int sliceSize = w * h;
//        int imageSize = sliceSize * d;
//        double data[] = image.getData();
//        double slice[] = new double[sliceSize];
//        ImageStack is = new ImageStack(w, h);
//
//        for (int i = 0; i < n; i++) {
//            int offset = i * imageSize;
//            for (int j = 0; j < d; j++) {
//                System.arraycopy(data, offset + j * sliceSize, slice, 0, sliceSize);
//
//                FloatProcessor processor = new FloatProcessor(w, h, slice);
//                is.addSlice(String.valueOf(i), processor);
//            }
//        }
//
//        imp.setStack(is);
    }

    public static ImageGeneric_ convertToXmipp(ImagePlus imp) {
        ImageGeneric_ image = new ImageGeneric_();

//        int w = imp.getWidth();
//        int h = imp.getHeight();
//        int d = imp.getStackSize();
//
//        double data[] = new double[w * h * d];
//        for (int i = 0; i < d; i++) {
//            float slice[] = (float[]) imp.getStack().getProcessor(i + 1).getPixels();
//            System.arraycopy(Tools.toDouble(slice), 0, data, i * w * h, w * h);
//        }
//        try {
//            image.setData(w, h, d, data);
//        } catch (Exception ex) {
//            ex.printStackTrace(System.err);
//            IJ.error(ex.getMessage());
//        }

        return image;
    }

    public static boolean saveImage(ImagePlus imp, String filename) {
//        try {
//            ImageGeneric image = convertToXmipp(imp);
//
//            image.write(filename);
//            return true;
//        } catch (Exception ex) {
//            ex.printStackTrace(System.err);
//            IJ.error(ex.getMessage());
//        }

        return false;
    }
}

abstract class ProcessorCreator_ {

    public abstract ImageProcessor getProcessor(ImageGeneric_ image, int slice) throws Exception;
}

class FloatProcessorCreator_ extends ProcessorCreator_ {

    public ImageProcessor getProcessor(ImageGeneric_ image, int slice) throws Exception {
        return new FloatProcessor(image.getXDim(), image.getYDim(), image.getArrayFloat(slice), null);
    }
}

class ShortProcessorCreator_ extends ProcessorCreator_ {

    public ImageProcessor getProcessor(ImageGeneric_ image, int slice) throws Exception {
        return new ShortProcessor(image.getXDim(), image.getYDim(), image.getArrayShort(slice), null);
    }
}

class ByteProcessorCreator_ extends ProcessorCreator_ {

    public ImageProcessor getProcessor(ImageGeneric_ image, int slice) throws Exception {
        return new ByteProcessor(image.getXDim(), image.getYDim(), image.getArrayByte(slice), null);
    }
}