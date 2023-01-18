#include <filesystem>
#include <iostream>
#include <vector>
namespace fs = std::filesystem;

int main()
{
    std::ios_base::sync_with_stdio(false);
    for (const auto &entry : fs::recursive_directory_iterator(".")) {
            std::cout << entry.path().string() << std::endl;
            
    }
    return 0;
}
