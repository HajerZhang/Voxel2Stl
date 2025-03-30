#include <voxel2stl.hpp>

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "用法: " << argv[0] << " <输入VTK文件> <输出STL文件> [细化指数] <等值面值> [平滑迭代次数]" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string inputVtk = argv[1];
    std::string outputStl = argv[2];
    int refineLevel = (argc > 3) ? std::stoi(argv[3]) : 1;
    if (refineLevel < 1) {
        std::cerr << "细化指数必须大于等于1！" << std::endl;
        return EXIT_FAILURE;
    }
    double isoValue = std::stod(argv[4]);
    int smoothIterations = (argc > 5) ? std::stoi(argv[4]) : 20;

    try {
        #ifdef IPDI
        Voxel2Stl *voxel2stl = new Voxel2Stl(inputVtk, refineLevel, isoValue, smoothIterations);
        voxel2stl->OutputStl(outputStl);
        delete voxel2stl;
        #endif
        #ifdef DERECT
        Voxel2Stl::VTK2Stl_D(inputVtk, outputStl, isoValue, smoothIterations);
        #endif
        #ifdef ONEGET
        Voxel2Stl::VTK2Stl_O(inputVtk, outputStl, isoValue, smoothIterations);
        #endif
    } catch (const std::exception& e) {
        std::cerr << "错误: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
