#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

int main(int argc, char *argv[])
{
	string str = "#include \"bitmorev5.h\"\n\nint main(int argc, char *argv[])\n{\n";
	/*
	constexpr size_t nitems = (1ULL << 14);
  constexpr size_t nbytes_per_row = 1024*15;
  constexpr size_t nservers = 5;//s+1
  constexpr size_t soundness = 128;
	*/
	str += "  constexpr size_t nitems = (1ULL <<" + to_string(14) + ");\n";
	str += "  constexpr size_t nbytes_per_row = "+to_string(1024*15)+";\n";
	str += "  constexpr size_t nservers = "+to_string(5)+";\n";
	str += "  constexpr size_t soundness = "+to_string(128)+";\n"
	ifstream inFile;
	inFile.open("runBitMoreBase");//open the input file

	stringstream strStream;
	strStream << inFile.rdbuf();//read the file
	str += strStream.str();//str holds the content of the file

	cout << str << endl;
	inFile.close();
	return 0;
}
