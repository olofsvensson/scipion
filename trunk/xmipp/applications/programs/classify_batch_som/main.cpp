/***************************************************************************
 *
 * Authors:     Alberto Pascual Montano (pascual@cnb.uam.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

// To avoid problems with long template names
#pragma warning(disable:4786)

#include <fstream>

#include <classification/batch_som.h>

/* Prototypes -============================================================= */

void Usage(char **argv);

/* Main function -============================================================= */

main(int argc, char** argv)
{


    /* Input Parameters ======================================================== */

    FileName       fn_in;     // input file
    FileName       fn_out;    // output file
    FileName       cb_in = "";    // Code vectors input file
    FileName       tmpN;  // Temporary variable
    unsigned       iter = 1000; // Iteration number
    unsigned       verb = 0; // Verbosity level
    bool           norm = 1; // Normalize?
    unsigned       xdim = 10; // X-dimension (-->)
    unsigned       ydim = 5; // Y-dimension
    float         radius_0; // Initial radius value
    string         layout = "HEXA"; // Type of layout (Topology)
    bool         saveClusters = false;    // Save clusters in separate files

    /* Parameters ============================================================== */
    try
    {

        fn_in = getParameter(argc, argv, "-i");

        if (checkParameter(argc, argv, "-o"))
            fn_out = getParameter(argc, argv, "-cout");
        else
        {
            Usage(argv);
            exit(EXIT_FAILURE);
        }

        if (checkParameter(argc, argv, "-cvin"))
            cb_in = getParameter(argc, argv, "-cvin");

        ydim = textToInteger(getParameter(argc, argv, "-ydim", "5"));
        xdim = textToInteger(getParameter(argc, argv, "-xdim", "10"));

        if (checkParameter(argc, argv, "-hexa"))
        {
            if (checkParameter(argc, argv, "-rect"))
            {
                cout << "Error: you can not define two topologies" << endl;
                exit(EXIT_FAILURE);
            }
            layout = "HEXA";
        }
        else if (checkParameter(argc, argv, "-rect"))
            layout = "RECT";


        if (checkParameter(argc, argv, "-radius"))
            radius_0 = textToFloat(getParameter(argc, argv, "-radius"));
        else
            if (xdim > ydim)
                radius_0 = xdim;
            else
                radius_0 = ydim;

        iter = textToInteger(getParameter(argc, argv, "-iter", "1000"));
        verb = textToInteger(getParameter(argc, argv, "-verb", "0"));

        if (checkParameter(argc, argv, "-norm"))
            norm = true;
        else norm = false;

        if (checkParameter(argc, argv, "-saveclusters"))
            saveClusters = true;
        else saveClusters = false;


        if (argc == 1)
        {
            Usage(argv);
        }

    }
    catch (Xmipp_error XE)
    {
        cout << XE;
        Usage(argv);
    }


    /* Some validations ===================================================== */


    if (iter < 1)
    {
        cerr << argv[0] << ": invalid value for iter (must be > 1): " << iter << endl;
        exit(EXIT_FAILURE);
    }

    if (verb < 0 || verb > 2)
    {
        cerr << argv[0] << ": invalid value for verbosity (must be between 0 and 2): " << verb << endl;
        exit(EXIT_FAILURE);
    }


    if (radius_0 < 1)
    {
        cerr << argv[0] << ": invalid value for radius (must be > 1): " << radius_0 << endl;
        exit(EXIT_FAILURE);
    }

    if (xdim < 1)
    {
        cerr << argv[0] << ": invalid value for xdim (must be > 1): " << xdim << endl;
        exit(EXIT_FAILURE);
    }

    if (ydim < 1)
    {
        cerr << argv[0] << ": invalid value for ydim (must be > 1): " << ydim << endl;
        exit(EXIT_FAILURE);
    }


    /* Shows parameters ===================================================== */

    cout << endl << "Parameters used: " << endl;
    cout << "Input data file : " << fn_in << endl;
    cout << "Output file name : " << fn_out << endl;
    if (cb_in != "")
        cout << "Input code vectors file name : " << cb_in << endl;
    if (saveClusters)
        cout << "Save clusters in separate files: " << fn_out << ".(cluster number)" << endl;
    cout << "Horizontal dimension (Xdim) = " << xdim << endl;
    cout << "Vertical dimension (Ydim) = " << ydim << endl;
    if (layout == "HEXA")
        cout << "Hexagonal topology " << endl;
    else
        cout << "Rectangular topology " << endl;

    cout << "Initial neighborhood radius (radius) = " << radius_0 << endl;

    cout << "Total number of iterations = " << iter << endl;
    cout << "verbosity level = " << verb << endl;
    if (norm)
        cout << "Normalize input data" << endl;
    else
        cout << "Do not normalize input data " << endl;


    /* Open training vector ================================================= */


    ifstream inStream(fn_in.c_str());
    if (!inStream)
    {
        cerr << argv[0] << ": can't open file " << fn_in << endl;
        exit(EXIT_FAILURE);
    }

    xmippCTVectors ts(0, true);

    cout << endl << "Reading data file " << fn_in << "....." << endl;
    inStream >> ts;


    /* Real stuff ============================================================== */


    try
    {
        if (norm)
        {
            cout << "Normalizing....." << endl;
            ts.normalize();        // Normalize input data
        }

        xmippMap* myMap;

        if (cb_in != "")
        {
            cout << "Reading codevectors file " << cb_in << "....." << endl;
            ifstream codeStream(cb_in.c_str());
            if (!codeStream)
            {
                cerr << argv[0] << ": can't open file " << cb_in << endl;
                exit(EXIT_FAILURE);
            }
            myMap = new xmippMap(codeStream);
        }
        else
            myMap = new xmippMap(layout, xdim, ydim, ts);

        xmippBatchSOM *thisSOM;
        Descent radius(radius_0, 1);            // radius decreases linearly to 1
        thisSOM = new xmippBatchSOM(radius, iter);    // Creates BatchSOM Algorithm

        xmippTextualListener myListener;     // Define the listener class
        myListener.setVerbosity() = verb;     // Set verbosity level
        thisSOM->setListener(&myListener);        // Set Listener

        thisSOM->train(*myMap, ts);                // Train algorithm

        // Test algorithm
        cout << endl;
        double dist = thisSOM->test(*myMap, ts);
        cout << endl << "Quantization error : " <<  dist << endl;

        // Classifying
        cout << "Classifying....." << endl;
        myMap->classify(&ts);

        // Calibrating
        cout << "Calibrating....." << endl;
        myMap->calibrate(ts);

        /*******************************************************
            Saving all kind of Information
        *******************************************************/

        cout << "Saving algorithm information as " << fn_out << ".inf ....." << endl;
        tmpN = fn_out.c_str() + (string) ".inf";
        ofstream infS(tmpN.c_str());
        infS << "Kohonen BatchSOM algorithm" << endl << endl;
        infS << "Input data file : " << fn_in << endl;
        if (cb_in != "")
            infS << "Input code vectors file : " << cb_in << endl;
        infS << "Code vectors output file : " << fn_out <<  ".cod" << endl;
        infS << "Algorithm information output file : " << fn_out <<  ".inf" << endl;
        infS << "Number of feature vectors: " << ts.size() << endl;
        infS << "Number of variables: " << ts.itemAt(0).size() << endl;
        infS << "Horizontal dimension (Xdim) = " << xdim << endl;
        infS << "Vertical dimension (Ydim) = " << ydim << endl;
        if (layout == "HEXA")
            infS << "Hexagonal topology " << endl;
        else
            infS << "Rectangular topology " << endl;
        infS << "Initial neighborhood radius (radius) = " << radius_0 << endl;
        infS << "Total number of iterations = " << iter << endl;
        if (norm)
            infS << "Input data normalized" << endl;
        else
            infS << "Input data not normalized" << endl;
        infS << "Quantization error : " <<  dist << endl;

        infS.flush();

        // assign data to neurons
        if (saveClusters)
        {
            cout << "Saving neurons assigments ....." << endl;
            for (unsigned i = 0; i < myMap->size(); i++)
            {
                tmpN = fn_out.c_str() + (string) "."  + ItoA(i);
                ofstream cStream(tmpN.c_str());
                for (int j = 0; j < myMap->classifAt(i).size(); j++)
                    cStream << myMap->classifAt(i)[j] << endl;
                cStream.flush();
            }
        }

        // save .vs file to be compatible with SOM_PAK
        cout << "Saving visual file as " << fn_out << ".vs ....." << endl;
        tmpN = fn_out.c_str() + (string) ".vs";
        ofstream vsStream(tmpN.c_str());
        vsStream << ts.theItems[0].size() << " " << myMap->layout() << " " << myMap->width() << " " << myMap->height() << " bubble" << endl;
        for (int i = 0; i < ts.size(); i++)
        {
            int j = myMap->winner(ts, i);
            vsStream << myMap->indexToPos(j).first << " " << myMap->indexToPos(j).second << " " << eDist(myMap->theItems[j], ts.theItems[i]) << " " << ts.theTargets[i] << endl;
        }
        vsStream.flush();

        // save .his file (Histogram)
        cout << "Saving code vectors histogram file as " << fn_out << ".his ....." << endl;
        tmpN = fn_out.c_str() + (string) ".his";
        ofstream hisStream(tmpN.c_str());
        myMap->printHistogram(hisStream);
        hisStream.flush();

        // save .err file (Average Quantization Error)
        cout << "Saving code vectors average quantization error file as " << fn_out << ".err ....." << endl;
        tmpN = fn_out.c_str() + (string) ".err";
        ofstream errStream(tmpN.c_str());
        myMap->printQuantError(errStream);
        errStream.flush();

        if (norm)
        {
            cout << "Denormalizing code vectors....." << endl;
            myMap->unNormalize(ts.getNormalizationInfo()); // de-normalize codevectors
        }

        cout << "Saving code vectors as " << fn_out << ".cod ....." << endl;
        tmpN = fn_out.c_str() + (string) ".cod";
        ofstream codS(tmpN.c_str());
        codS << *myMap;
        codS.flush();

        cout << endl;

        delete myMap;

    }
    catch (const exception& e)
    {
        cout << e.what() << endl;
    }
    return 0;
}


/* ------------------------------------------------------------------------- */
/* Help Message for this Program                                             */
/* ------------------------------------------------------------------------- */
void Usage(char **argv)
{
    printf(
        "\nUsage: %s [Purpose and Parameters]"
        "\nPurpose: Kohonen Self-Organizing Feature Maps (Batch Training)"
        "\n"
        "\nParameter Values: (note space before value)"
        "\n"
        "\n    -i      file_in           Input data file (plain data)"
        "\n    -o      file_out          Base name for output data files:"
        "\n    -cvin    file_in          Codevectors input file"
        "\n    -saveclusters           save clusters in separate files (Default = No)"
        "\n    -xdim   H-dimension       Horizontal size of the map (default = 10)"
        "\n    -ydim   V-dimension       Vertical size of the map (default = 5)"
        "\n    -hexa               Hexagonal topology (default)"
        "\n    -rect               Rectangular topology"
        "\n    -radius    radius        Initial neighborhood radius"
        "\n             (default = max(xdim, ydim)"
        "\n    -iter      iterations     Number of iterations (default = 1000)"
        "\n    -norm                   Normalize training data (default)"
        "\n    -verb      verbosity      Information level while running: "
        "\n             0: No information (default)"
        "\n             1: Progress bar"
        "\n             2: Code vectors change between iterations"
        "\n      \n"
        , argv[0]);
    exit(0);
}

