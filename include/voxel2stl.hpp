#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkContourFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSTLWriter.h>
#include <vtkCellDataToPointData.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkFloatArray.h>
#include <vtkFeatureEdges.h>
#include <vtkAppendPolyData.h>
#include <vtkStructuredPointsWriter.h>
#include <iostream>
#include <stdexcept>

#define MODE 1
#if MODE == 1
#define IPDI
#elif MODE == 2
#define DERECT
#elif MODE == 3
#define ONEGET
#endif


class Voxel2Stl
{
    private:

    int m_dims_voxel[3];
    double m_spacing_voxel[3];
    double m_origin_voxel[3];

    int m_numCells;
    std::vector<double> m_cellScalars;

    int m_refineLevel;

    int m_numPoints;
    std::vector<double> m_pointScalars;

    int m_smoothIterations;
    double m_isoValue;

    void CreateRefinedPoints();

    public:
    Voxel2Stl();
    Voxel2Stl(const std::string& inputVtkFile, int refineLevel, double isoValue, int smoothIterations = 20);
    ~Voxel2Stl();

    void OutputStl(const std::string& outputStlFile);
    static void VTK2Stl_D(const std::string& inputVtkFile, const std::string& outputStlFile, double isoValue, int smoothIterations = 20);
    static void VTK2Stl_O(const std::string& inputVtkFile, const std::string& outputStlFile, double isoValue, int smoothIterations = 20);
};

