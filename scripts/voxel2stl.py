import vtk
import os
import numpy as np

def vtk_to_smooth_stl(input_vtk_file, output_stl_file, iso_value, smooth_iterations=20):
    if not os.path.exists(input_vtk_file):
        raise FileNotFoundError(f"文件 {input_vtk_file} 不存在，请检查路径！")
    
    reader = vtk.vtkStructuredPointsReader()
    reader.SetFileName(input_vtk_file)
    reader.Update()
    
    output = reader.GetOutput()
    if output is None or output.GetNumberOfCells() == 0:
        raise ValueError("输入 VTK 文件无有效数据或未正确解析！")
    
    scalars = output.GetCellData().GetScalars()
    if scalars is None:
        raise ValueError("VTK 文件中未找到标量数据！")
    
    dims = output.GetDimensions()
    print(f"原始模型维度：{dims}")
    
    num_cells = scalars.GetNumberOfTuples()
    scalar_data = np.zeros(num_cells)
    for i in range(num_cells):
        scalar_data[i] = scalars.GetValue(i)
    
    scalar_data = scalar_data.reshape(dims[2] - 1, dims[1] - 1, dims[0] - 1)
    
    padded_data = np.pad(scalar_data, pad_width=1, mode='constant', constant_values=0)
    print(f"扩展后的模型维度：{padded_data.shape}")
    
    new_dims = (padded_data.shape[2] + 1, padded_data.shape[1] + 1, padded_data.shape[0] + 1)
    new_structured_points = vtk.vtkStructuredPoints()
    new_structured_points.SetDimensions(new_dims)
    new_structured_points.SetOrigin(output.GetOrigin()-np.array([1, 1, 1]) * output.GetSpacing())
    new_structured_points.SetSpacing(output.GetSpacing())
    
    new_scalars = vtk.vtkDoubleArray()
    new_scalars.SetNumberOfValues(padded_data.size)
    new_scalars.SetName("PaddedScalars")
    flat_data = padded_data.ravel()
    for i in range(flat_data.size):
        new_scalars.SetValue(i, flat_data[i])
    
    new_structured_points.GetCellData().SetScalars(new_scalars)
    
    cell_to_point = vtk.vtkCellDataToPointData()
    cell_to_point.SetInputData(new_structured_points)
    cell_to_point.Update()
    output_converted = cell_to_point.GetOutput()
    
    scalars = output_converted.GetPointData().GetScalars()
    if scalars is None:
        raise ValueError("转换后仍未找到标量数据！")
    
    min_val, max_val = scalars.GetRange()
    print(f"标量值范围：最小值={min_val}, 最大值={max_val}")
    
    if not (min_val <= iso_value <= max_val):
        raise ValueError(f"指定的等值面值 {iso_value} 不在标量值范围内！")
    
    contour_filter = vtk.vtkContourFilter()
    contour_filter.SetInputData(output_converted)
    contour_filter.SetValue(0, iso_value)
    contour_filter.Update()
    
    if contour_filter.GetOutput().GetNumberOfCells() == 0:
        raise ValueError("等值面提取失败，可能是等值面值不合适！")
    
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputConnection(contour_filter.GetOutputPort())
    smoother.SetNumberOfIterations(smooth_iterations)
    smoother.SetRelaxationFactor(0.1)  # 控制平滑程度，默认0.1
    smoother.FeatureEdgeSmoothingOff()  # 可根据需要启用边界平滑
    smoother.BoundarySmoothingOn()
    smoother.Update()
    
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetInputConnection(smoother.GetOutputPort())
    stl_writer.SetFileName(output_stl_file)
    stl_writer.Write()
    
    print(f"平滑后的等值面 {iso_value} 已保存到 {output_stl_file} 中。")

input_vtk = "Input.vtk" 
output_stl = "Output.stl" 
iso_value = 0.9
smooth_iterations = 0

vtk_to_smooth_stl(input_vtk, output_stl, iso_value, smooth_iterations)