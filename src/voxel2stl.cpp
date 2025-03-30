#include <voxel2stl.hpp>

Voxel2Stl::Voxel2Stl() {
    
}

Voxel2Stl::~Voxel2Stl() {
    m_cellScalars.clear();
    m_pointScalars.clear();
}

Voxel2Stl::Voxel2Stl(const std::string& inputVtkFile, int refineLevel, double isoValue, int smoothIterations) {
    
    auto reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader->SetFileName(inputVtkFile.c_str());
    reader->Update();
    
    vtkStructuredPoints* output = reader->GetOutput();
    if (!output || output->GetNumberOfCells() == 0) {
        throw std::runtime_error("输入 VTK 文件无有效数据或未正确解析！");
    }

    vtkDataArray* cellScalars = output->GetCellData()->GetScalars();
    if (!cellScalars) {
        throw std::runtime_error("未找到 CellData 的标量数据！");
    }

    double minVal = cellScalars->GetRange()[0];
    double maxVal = cellScalars->GetRange()[1];
    std::cout << "标量值范围：最小值=" << minVal << ", 最大值=" << maxVal << std::endl;
    if (!(minVal <= isoValue && isoValue <= maxVal)) {
        throw std::runtime_error("指定的等值面值不在标量值范围内！");
    }

    int *dims = output->GetDimensions();
    m_dims_voxel[0] = dims[0] - 1;
    m_dims_voxel[1] = dims[1] - 1;
    m_dims_voxel[2] = dims[2] - 1;

    double *spacing = output->GetSpacing();
    m_spacing_voxel[0] = spacing[0];
    m_spacing_voxel[1] = spacing[1];
    m_spacing_voxel[2] = spacing[2];

    double *origin = output->GetOrigin();
    m_origin_voxel[0] = origin[0];
    m_origin_voxel[1] = origin[1];
    m_origin_voxel[2] = origin[2];

    m_numCells = output->GetNumberOfCells();
    m_cellScalars.resize(m_numCells);
    for (int i = 0; i < m_numCells; i++) {
        m_cellScalars[i] = cellScalars->GetComponent(i, 0);
    }

    m_refineLevel = refineLevel;
    m_isoValue = isoValue;
    m_smoothIterations = smoothIterations;
}

void Voxel2Stl::CreateRefinedPoints()
{
    int refinedDims[3];
    refinedDims[0] = m_dims_voxel[0] * m_refineLevel;
    refinedDims[1] = m_dims_voxel[1] * m_refineLevel;
    refinedDims[2] = m_dims_voxel[2] * m_refineLevel;

    m_numPoints = (refinedDims[0] + 1) * (refinedDims[1] + 1) * (refinedDims[2] + 1);
    m_pointScalars.resize(m_numPoints);

    double refinedSpacing[3];
    refinedSpacing[0] = m_spacing_voxel[0] / m_refineLevel;
    refinedSpacing[1] = m_spacing_voxel[1] / m_refineLevel;
    refinedSpacing[2] = m_spacing_voxel[2] / m_refineLevel;

    double radius = 1.732;
    double sigma = radius / 3.0;

    for(int k = 0; k < refinedDims[2] + 1; k++) {
        for(int j = 0; j < refinedDims[1] + 1; j++) {
            for(int i = 0; i < refinedDims[0] + 1; i++) {

                int idx = k * (refinedDims[0] + 1) * (refinedDims[1] + 1) + j * (refinedDims[0] + 1) + i;

                double x_point_grid = i / (1.0 * m_refineLevel);
                double y_point_grid = j / (1.0 * m_refineLevel);
                double z_point_grid = k / (1.0 * m_refineLevel);

                int neighbor_cell_xMin = std::max(0, (int)std::ceil(x_point_grid - radius));
                int neighbor_cell_xMax = std::min(m_dims_voxel[0] - 1, (int)std::floor(x_point_grid + radius));
                int neighbor_cell_yMin = std::max(0, (int)std::ceil(y_point_grid - radius));
                int neighbor_cell_yMax = std::min(m_dims_voxel[1] - 1, (int)std::floor(y_point_grid + radius));
                int neighbor_cell_zMin = std::max(0, (int)std::ceil(z_point_grid - radius));
                int neighbor_cell_zMax = std::min(m_dims_voxel[2] - 1, (int)std::floor(z_point_grid + radius));

                double value = 0.0;
                double weightSum = 0.0;
                for(int k2 = neighbor_cell_zMin; k2 <= neighbor_cell_zMax; k2++) {
                    for(int j2 = neighbor_cell_yMin; j2 <= neighbor_cell_yMax; j2++) {
                        for(int i2 = neighbor_cell_xMin; i2 <= neighbor_cell_xMax; i2++) {

                            int idx_neighbor = k2 * m_dims_voxel[0] * m_dims_voxel[1] + j2 * m_dims_voxel[0] + i2;

                            double x_center = (i2 + 0.5);
                            double y_center = (j2 + 0.5);
                            double z_center = (k2 + 0.5);

                            double distance = sqrt(pow(x_point_grid - x_center, 2) + pow(y_point_grid - y_center, 2) + pow(z_point_grid - z_center, 2));
                            #define GAUSS
                            #ifdef GAUSS
                            double factor = exp(-distance * distance / (sigma * sigma)); 
                            #else
                            double factor = radius - distance;
                            #endif

                            if(factor > 0) {
                                #ifdef GAUSS
                                double weight = factor;
                                #else
                                double weight = factor * factor * factor;
                                #endif
                                value += weight * m_cellScalars[idx_neighbor];
                                weightSum += weight;
                            }
                        }
                    }
                }

                if(weightSum > 0) {
                    value /= weightSum;
                } else {
                    value = 0.0;
                }
                m_pointScalars[idx] = value;
            }
        }
    }

}

void Voxel2Stl::OutputStl(const std::string& outputStlFile) 
{
    CreateRefinedPoints();
    int refinedDims[3];
    refinedDims[0] = m_dims_voxel[0] * m_refineLevel + 1;
    refinedDims[1] = m_dims_voxel[1] * m_refineLevel + 1;
    refinedDims[2] = m_dims_voxel[2] * m_refineLevel + 1;
    double refinedSpacing[3];
    refinedSpacing[0] = m_spacing_voxel[0] / m_refineLevel;
    refinedSpacing[1] = m_spacing_voxel[1] / m_refineLevel;
    refinedSpacing[2] = m_spacing_voxel[2] / m_refineLevel;

    int paddedDims[3];
    paddedDims[0] = refinedDims[0] + 2;
    paddedDims[1] = refinedDims[1] + 2;
    paddedDims[2] = refinedDims[2] + 2;

    // Padded Points
    auto refinedPoints = vtkSmartPointer<vtkStructuredPoints>::New();
    refinedPoints->SetDimensions(paddedDims[0], paddedDims[1], paddedDims[2]);
    refinedPoints->SetSpacing(m_spacing_voxel[0] / m_refineLevel, m_spacing_voxel[1] / m_refineLevel, m_spacing_voxel[2] / m_refineLevel);
    refinedPoints->SetOrigin(m_origin_voxel[0] - m_spacing_voxel[0] / m_refineLevel, m_origin_voxel[1] - m_spacing_voxel[1] / m_refineLevel, m_origin_voxel[2] - m_spacing_voxel[2] / m_refineLevel);

    auto refinedScalars = vtkSmartPointer<vtkDoubleArray>::New();
    refinedScalars->SetNumberOfValues(paddedDims[0] * paddedDims[1] * paddedDims[2]);
    refinedScalars->FillComponent(0, 0);
    refinedScalars->SetName("PaddedScalars");

    std::vector<double> paddingScalars(paddedDims[0] * paddedDims[1] * paddedDims[2], 0.0);

    for(int k = 0; k < refinedDims[2]; k++) {
        for(int j = 0; j < refinedDims[1]; j++) {
            for(int i = 0; i < refinedDims[0]; i++) {
                int idx_Origin = k * refinedDims[0] * refinedDims[1] + j * refinedDims[0] + i;
                int idx_Padded = (k + 1) * paddedDims[0] * paddedDims[1] + (j + 1) * paddedDims[0] + (i + 1);
                paddingScalars[idx_Padded] = m_pointScalars[idx_Origin];
            }
        }
    }
    
    int size = paddedDims[0] * paddedDims[1] * paddedDims[2];
    for(int i = 0; i < size; i++){
        refinedScalars->SetTuple1(i, paddingScalars[i]);
    }

    refinedPoints->GetPointData()->SetScalars(refinedScalars);
    refinedPoints->GetPointData()->SetActiveScalars("PaddedScalars");
    refinedPoints->GetPointData()->SetActiveAttribute("PaddedScalars", vtkDataSetAttributes::SCALARS);

    // // Save the refined points to a VTK file for debugging
    // vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    // writer->SetFileName("output2.vtk");
    // writer->SetInputData(refinedPoints);
    // writer->Write();

    auto contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputData(refinedPoints);
    contourFilter->SetValue(0, m_isoValue);
    contourFilter->Update();

    if (contourFilter->GetOutput()->GetNumberOfCells() == 0) {
        throw std::runtime_error("等值面提取失败，可能是等值面值不合适！");
    }

    auto smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(contourFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(m_smoothIterations);
    smoothFilter->SetRelaxationFactor(0.5);
    smoothFilter->FeatureEdgeSmoothingOn();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    auto stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetInputConnection(smoothFilter->GetOutputPort());
    stlWriter->SetFileName(outputStlFile.c_str());
    stlWriter->SetFileTypeToBinary();
    stlWriter->Write();

    if(stlWriter->GetErrorCode() != 0) {
        throw std::runtime_error("STL 文件写入失败！");
    }else{
        std::cout << "STL 文件写入成功 "<< outputStlFile << std::endl;
    }
}

void Voxel2Stl::VTK2Stl_D(const std::string& inputVtkFile, const std::string& outputStlFile, double isoValue, int smoothIterations) 
{
    // 读取 VTK 文件
    auto reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader->SetFileName(inputVtkFile.c_str());
    reader->Update();
    
    vtkStructuredPoints* output = reader->GetOutput();
    if (!output || output->GetNumberOfCells() == 0) {
        throw std::runtime_error("输入 VTK 文件无有效数据或未正确解析！");
    }

    // 获取标量数据范围
    vtkDataArray* scalars = output->GetCellData()->GetScalars();
    if (!scalars) {
        throw std::runtime_error("转换后仍未找到标量数据！");
    }
    
    double minVal = scalars->GetRange()[0];
    double maxVal = scalars->GetRange()[1];
    std::cout << "标量值范围：最小值=" << minVal << ", 最大值=" << maxVal << std::endl;
    
    if (!(minVal <= isoValue && isoValue <= maxVal)) {
        throw std::runtime_error("指定的等值面值不在标量值范围内！");
    }

    auto PaddledCell = vtkSmartPointer<vtkStructuredPoints>::New();
    PaddledCell->SetDimensions(output->GetDimensions()[0] + 2, output->GetDimensions()[1] + 2, output->GetDimensions()[2] + 2);
    PaddledCell->SetSpacing(output->GetSpacing()[0], output->GetSpacing()[1], output->GetSpacing()[2]);
    PaddledCell->SetOrigin(output->GetOrigin()[0] - output->GetSpacing()[0], output->GetOrigin()[1] - output->GetSpacing()[1], output->GetOrigin()[2] - output->GetSpacing()[2]);

    auto PaddledScalars = vtkSmartPointer<vtkDoubleArray>::New();
    PaddledScalars->SetNumberOfValues((output->GetDimensions()[0] + 2) * (output->GetDimensions()[1] + 2) * (output->GetDimensions()[2] + 2));
    PaddledScalars->FillComponent(0, 0);
    PaddledScalars->SetName("PaddledScalars");

    std::vector<double> paddingScalars((output->GetDimensions()[0] + 2) * (output->GetDimensions()[1] + 2) * (output->GetDimensions()[2] + 2), 0.0);
    int dims[3];
    dims[0] = output->GetDimensions()[0] - 1;
    dims[1] = output->GetDimensions()[1] - 1;
    dims[2] = output->GetDimensions()[2] - 1;

    for(int k = 0; k < dims[2]; k++) {
        for(int j = 0; j < dims[1]; j++) {
            for(int i = 0; i < dims[0]; i++) {
                int idx_Origin = k * dims[0] * dims[1] + j * dims[0] + i;
                int idx_Padded = (k + 1) * (dims[0] + 2) * (dims[1] + 2) + (j + 1) * (dims[0] + 2) + (i + 1);
                paddingScalars[idx_Padded] = scalars->GetComponent(idx_Origin, 0);
            }
        }
    }

    for(int i = 0; i < (output->GetDimensions()[0] + 1) * (output->GetDimensions()[1] + 1) * (output->GetDimensions()[2] + 1); i++){
        PaddledScalars->SetTuple1(i, paddingScalars[i]);
    }

    PaddledCell->GetCellData()->SetScalars(PaddledScalars);
    PaddledCell->GetCellData()->SetActiveScalars("PaddledScalars");
    PaddledCell->GetCellData()->SetActiveAttribute("PaddledScalars", vtkDataSetAttributes::SCALARS);
    
    // CELL_DATA 转换为 POINT_DATA
    auto cellToPoint = vtkSmartPointer<vtkCellDataToPointData>::New();
    cellToPoint->SetInputData(PaddledCell);
    cellToPoint->Update();

    vtkImageData* outputConverted = vtkImageData::SafeDownCast(cellToPoint->GetOutput());
    if (!outputConverted) {
        throw std::runtime_error("转换后的数据无法转换为 vtkImageData！");
    }

    // 提取等值面
    auto contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputData(outputConverted);
    contourFilter->SetValue(0, isoValue);
    contourFilter->Update();
    
    if (contourFilter->GetOutput()->GetNumberOfCells() == 0) {
        throw std::runtime_error("等值面提取失败，可能是等值面值不合适！");
    }
    
    // 平滑处理
    auto smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoother->SetInputConnection(contourFilter->GetOutputPort());
    smoother->SetNumberOfIterations(smoothIterations);
    smoother->SetRelaxationFactor(0.5);
    smoother->FeatureEdgeSmoothingOff();
    smoother->BoundarySmoothingOn();
    smoother->Update();
    
    // 写入 STL 文件
    auto stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetInputConnection(smoother->GetOutputPort());
    stlWriter->SetFileName(outputStlFile.c_str());
    stlWriter->Write();
    
    std::cout << "平滑后的等值面 " << isoValue << " 已保存到 " << outputStlFile << " 中。" << std::endl;
}

void Voxel2Stl::VTK2Stl_O(const std::string& inputVtkFile, const std::string& outputStlFile, double isoValue, int smoothIterations) 
{
    // 读取 VTK 文件
    auto reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    reader->SetFileName(inputVtkFile.c_str());
    reader->Update();
    
    vtkStructuredPoints* output = reader->GetOutput();
    if (!output || output->GetNumberOfCells() == 0) {
        throw std::runtime_error("输入 VTK 文件无有效数据或未正确解析！");
    }

    // 获取标量数据范围
    vtkDataArray* scalars = output->GetCellData()->GetScalars();
    if (!scalars) {
        throw std::runtime_error("转换后仍未找到标量数据！");
    }
    
    double minVal = scalars->GetRange()[0];
    double maxVal = scalars->GetRange()[1];
    std::cout << "标量值范围：最小值=" << minVal << ", 最大值=" << maxVal << std::endl;
    
    if (!(minVal <= isoValue && isoValue <= maxVal)) {
        throw std::runtime_error("指定的等值面值不在标量值范围内！");
    }
    

    int dims[3];
    dims[0] = output->GetDimensions()[0];
    dims[1] = output->GetDimensions()[1];
    dims[2] = output->GetDimensions()[2];

    std::vector<double> PointScalers(dims[0] * dims[1] * dims[2], 0.0);
    for(int k = 0; k < dims[2]; k++) {
        for(int j = 0; j < dims[1]; j++) {
            for(int i = 0; i < dims[0]; i++) {
                int idx = k * dims[0] * dims[1] + j * dims[0] + i;
                
                for(int l = -1; l < 1; l++) {
                    for(int m = -1; m < 1; m++) {
                        for(int n = -1; n < 1; n++) {
                            int Xcell = i + l;
                            int Ycell = j + m;
                            int Zcell = k + n;
                            if(Xcell >= 0 && Xcell < (dims[0] - 1) && Ycell >= 0 && Ycell < (dims[1] - 1) && Zcell >= 0 && Zcell < (dims[2] - 1)) {
                                int idx_neighbor = Zcell * (dims[0] - 1) * (dims[1] - 1) + Ycell * (dims[0] - 1) + Xcell;
                                double value = scalars->GetComponent(idx_neighbor, 0);
                                if(value > isoValue) {
                                    PointScalers[idx] = 1.0;
                                } else {
                                    PointScalers[idx] = 0.0;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // auto OnePoints = vtkSmartPointer<vtkStructuredPoints>::New();
    // OnePoints->SetDimensions(output->GetDimensions()[0], output->GetDimensions()[1], output->GetDimensions()[2]);
    // OnePoints->SetSpacing(output->GetSpacing()[0], output->GetSpacing()[1], output->GetSpacing()[2]);
    // OnePoints->SetOrigin(output->GetOrigin()[0], output->GetOrigin()[1], output->GetOrigin()[2]);

    // auto OneScalars = vtkSmartPointer<vtkDoubleArray>::New();
    // OneScalars->SetNumberOfValues(output->GetDimensions()[0] * output->GetDimensions()[1] * output->GetDimensions()[2]);
    // OneScalars->FillComponent(0, 0);
    // OneScalars->SetName("OneScalars");

    // int size = output->GetDimensions()[0] * output->GetDimensions()[1] * output->GetDimensions()[2];
    // for(int i = 0; i < size; i++){
    //     OneScalars->SetTuple1(i, PointScalers[i]);
    // }

    // OnePoints->GetPointData()->SetScalars(OneScalars);
    // OnePoints->GetPointData()->SetActiveScalars("OneScalars");
    // OnePoints->GetPointData()->SetActiveAttribute("OneScalars", vtkDataSetAttributes::SCALARS);

    // // Save the refined points to a VTK file for debugging
    // vtkSmartPointer<vtkStructuredPointsWriter> writer = vtkSmartPointer<vtkStructuredPointsWriter>::New();
    // writer->SetFileName("output2.vtk");
    // writer->SetInputData(OnePoints);
    // writer->Write();

    auto PaddledCell = vtkSmartPointer<vtkStructuredPoints>::New();
    PaddledCell->SetDimensions(output->GetDimensions()[0] + 2, output->GetDimensions()[1] + 2, output->GetDimensions()[2] + 2);
    PaddledCell->SetSpacing(output->GetSpacing()[0], output->GetSpacing()[1], output->GetSpacing()[2]);
    PaddledCell->SetOrigin(output->GetOrigin()[0] - output->GetSpacing()[0], output->GetOrigin()[1] - output->GetSpacing()[1], output->GetOrigin()[2] - output->GetSpacing()[2]);

    auto PaddledScalars = vtkSmartPointer<vtkDoubleArray>::New();
    PaddledScalars->SetNumberOfValues((output->GetDimensions()[0] + 2) * (output->GetDimensions()[1] + 2) * (output->GetDimensions()[2] + 2));
    PaddledScalars->FillComponent(0, 0);
    PaddledScalars->SetName("PaddledScalars");

    for(int k = 0; k < output->GetDimensions()[2]; k++) {
        for(int j = 0; j < output->GetDimensions()[1]; j++) {
            for(int i = 0; i < output->GetDimensions()[0]; i++) {
                int idx_Origin = k * output->GetDimensions()[0] * output->GetDimensions()[1] + j * output->GetDimensions()[0] + i;
                int idx_Padded = (k + 1) * (output->GetDimensions()[0] + 2) * (output->GetDimensions()[1] + 2) + (j + 1) * (output->GetDimensions()[0] + 2) + (i + 1);
                PaddledScalars->SetTuple1(idx_Padded, PointScalers[idx_Origin]);
            }
        }
    }

    PaddledCell->GetPointData()->SetScalars(PaddledScalars);
    PaddledCell->GetPointData()->SetActiveScalars("PaddledScalars");
    PaddledCell->GetPointData()->SetActiveAttribute("PaddledScalars", vtkDataSetAttributes::SCALARS);

    auto contourFilter = vtkSmartPointer<vtkContourFilter>::New();
    contourFilter->SetInputData(PaddledCell);
    contourFilter->SetValue(0, 1);
    contourFilter->Update();

    if (contourFilter->GetOutput()->GetNumberOfCells() == 0) {
        throw std::runtime_error("等值面提取失败，可能是等值面值不合适！");
    }

    auto smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(contourFilter->GetOutputPort());
    smoothFilter->SetNumberOfIterations(smoothIterations);
    smoothFilter->SetRelaxationFactor(0.5);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    auto stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
    stlWriter->SetInputConnection(smoothFilter->GetOutputPort());
    stlWriter->SetFileName(outputStlFile.c_str());
    stlWriter->SetFileTypeToBinary();
    stlWriter->Write();

    if(stlWriter->GetErrorCode() != 0) {
        throw std::runtime_error("STL 文件写入失败！");
    }else{
        std::cout << "STL 文件写入成功 "<< outputStlFile << std::endl;
    }

}





