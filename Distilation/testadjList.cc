#include <filesystem>

namespace recursive_directory_iterator = std::filesystem::recursive_directory_iterator;

int main (int argc, char *argv[])
{
for (auto& dirEntry : recursive_directory_iterator("../instances/mc/")){

}
  
  return 0;
}

