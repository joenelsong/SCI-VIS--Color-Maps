#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>


// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X,
                              const float *Y, const float *F)
{
        float value = 0;
        int cellLocation[2] = {0, 0};

        if ((pt[0] < X[dims[0]-1]) && (pt[0] >= X[0]) && (pt[1] < Y[dims[1]-1]) && (pt[1] >= Y[0]))
        {
                // Get x-coord for containing cell
                for (int i = 0; i < dims[0]-1; i++)
                {
                        if ((pt[0] > X[i]) && (pt[0] < X[i+1]))
                        {
                                cellLocation[0] = i;
                                break;
                        }
                }

                // Get y-coord for containing cell
                for (int j = 0; j < dims[1]-1; j++)
                {
                        if ((pt[1] > Y[j]) && (pt[1] < Y[j+1]))
                        {
                                cellLocation[1] = j;
                                break;
                        }
                }

                // Store scalars for corners of cell
                float FBbox[4];
                FBbox[0] = F[GetPointIndex(cellLocation, dims)];
                FBbox[1] = F[GetPointIndex(cellLocation, dims)+1];
                FBbox[2] = F[GetPointIndex(cellLocation, dims)+dims[0]];
                FBbox[3] = F[GetPointIndex(cellLocation, dims)+dims[0]+1];

                // Perform bilinear interpolation
                float valA = FBbox[0] + ((pt[0]-X[cellLocation[0]]) / (X[cellLocation[0]+1]-X[cellLocation[0]])) * (FBbox[1]-FBbox[0]);
                float valB = FBbox[2] + ((pt[0]-X[cellLocation[0]]) / (X[cellLocation[0]+1]-X[cellLocation[0]])) * (FBbox[3]-FBbox[2]);
                value = valA + ((pt[1]-Y[cellLocation[1]]) / (Y[cellLocation[1]+1]-Y[cellLocation[1]]))  * (valB-valA);
        }
        return value;
}

/*
float LinearInterpolate(float FofA, float FofB, float t)
{
    float ret;
    ret = FofA + t*(FofB-FofA);
    return ret;
}


float
EvaluateFieldAtLocation(const float *pt, const int *dims, 
                        const float *X, const float *Y, const float *F)
{
    // General Interpolation Equation: 
    // t = (X-A) / (B-A)
    // F(X) = F(A) + t*(F(B)-F(A))


    //printf("FLOAT %f, %f\n ", pt[0], pt[1]);
    int tempval = 0;

    float xaboveinterp;
    float xbelowinterp;
    float retinterp;

    // Interpolate Xabove point
    int i = 0;
    int yidxabove=-99;
    int yidxbelow=-99;
    for (i = 0; i < dims[1]; i++)
    {
        
        if (Y[i] >= pt[1]) {
            //printf("Y[i] >= pt[1] \n");
            yidxabove = i;
            yidxbelow = i-1;
            break;
        }
    }
    
    int xidxabove =-99;
    int xidxbelow = -99;
    for (i = 0; i < dims[0]; i++)
    {
        if (X[i] >= pt[0]) {
            xidxabove = i;
            xidxbelow = i-1;
            break;
        }
    }

    // Out of Bounds Check, return 0 if cant interpolate because of out of bounds
    if (yidxabove == -99 || yidxbelow == -99 || xidxbelow == -99 || xidxabove == -99) {
        printf("should be 0");
	return 0.0;
    }

    int FofAidx[2];
    FofAidx[0] = xidxabove;
    FofAidx[1] = yidxbelow;

    int FofBidx[2];
    FofBidx[0] = xidxabove;
    FofBidx[1] = yidxabove;


    printf("Dims: %d X %d \n", dims[0], dims[1]);
    printf("Z = %.2f X %.2f \n", pt[0], pt[1]);
    printf("Y[idxbelow] = Y[%d,%d] = %.3f,.%.3f \n", yidxbelow, yidxbelow, X[xidxbelow], Y[yidxbelow]);
    printf("Y[idxabove] = Y[%d,%d] = %.3f,.%.3f \n", yidxabove, yidxabove, X[xidxabove], Y[yidxabove]);


    float t = (pt[1]-Y[yidxbelow]) / (Y[yidxabove]-Y[yidxbelow]);

    xaboveinterp = LinearInterpolate( F[GetPointIndex(FofAidx, dims)], F[GetPointIndex(FofBidx, dims)], t );
    //printf("xaboveinterp Interpolated value: %f\n", xaboveinterp);

    // Interpolate XBelow point

    FofAidx[0] = xidxbelow;
    FofAidx[1] = yidxbelow;

    FofBidx[0] = xidxabove;
    FofBidx[1] = yidxbelow;

    xbelowinterp = LinearInterpolate( F[GetPointIndex(FofAidx, dims)], F[GetPointIndex(FofBidx, dims)], t );
    //printf("xbelowinterp Interpolated value: %f\n", xbelowinterp);


    // Interpolate X (finally)

    t = (pt[0] - X[xidxbelow]) / (X[xidxabove] - X[xidxbelow]);
    
    retinterp = LinearInterpolate( xbelowinterp, xaboveinterp, t );

    return retinterp;
}
*/

void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB)
{
// F has been normalized in the main    
    RGB[0] = F*255;
    RGB[1] = F*255;
    RGB[2] = 128 + (F*127);

}


// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB)
{
    if (F < 0.5) {
        RGB[0] = 2*F*255;
        RGB[1] = 2*F*255;
        RGB[2] = 128 + (2*F*127);
    }
    else {
        RGB[0] = 255 - (F-0.5)*255;
        RGB[1] = 255 - 2*(F-0.5)*255;
        RGB[2] = 255 - 2*(F-0.5)*255;
    }
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB)
{
    float s = 1.0;
    float v = 1.0;
    float h = F*360/60;

    //int cc = 6; // number of colors, 6 seems to be the most common.
    int i = (int)h;
    float f = h - i;
    float p = v * (1 - s);
    float q = v * (1 - f * s);
    float t = v * (1 - (1 - f) * s); 
    
    float r, g, b;
    
    switch(i) {
	case 0: r = v; g = t; b = p; break;
	case 1: r = q; g = v; b = p; break;
	case 2: r = p; g = v; b = t; break;
	case 3: r = p; g = q; b = v; break;
	case 4: r = t; g = p; b = v; break;
	case 5: r = v; g = p; b = q; break;
    }
    RGB[0] = r*255;
    RGB[1] = g*255;
    RGB[2] = b*255;
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
 //printf("%f  ",X[0]); printf("%f   ",X[9]); printf("%f   ",X[19]); printf("%f   ",X[29]); printf("%f   ",X[39]); printf("%f \n",X[49]);
 //printf("%d x %d x %d\n",dims[0], dims[1], dims[2]);
    
    int nx = 500;
    int ny = 500;
    float minF = 1.2; //jn
    float maxF = 5.02; //jn

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

    for (i = 0 ; i < nx ; i++)
        for (j = 0 ; j < ny ; j++)
        {
            // ITERATE OVER PIXELS
            float pt[2];
            
            pt[0] = -9 + (i*18) / (float)(nx-1);
            pt[1] = -9 + (j*18) / (float)(ny-1);
	    //printf("%f, %f \n", pt[0], pt[1]);
            float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);

            float normalizedF = (f-minF) / (maxF-minF); //...; see step 5 re 1.2->5.02
            
            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
            ApplyDifferenceColorMap(normalizedF, buffer[1]+offset);
            ApplyHSVColorMap(normalizedF, buffer[2]+offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
	cout << "test" << endl;
}
