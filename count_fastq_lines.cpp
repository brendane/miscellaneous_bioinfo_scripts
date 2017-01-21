#include<fstream>
#include<iostream>
#include<string>

void count_lines(std::istream& file) {
    std::string line;
    unsigned long count = 0;
    int i = 0;
    while(std::getline(file, line)) {
        if(line[0] == '@' && i == 0) {
            count += 1;
        } else if(i == 0) {
            std::cerr << line << std::endl;
            std::cerr << "Incorrect formatting" << std::endl;
            throw 1;
        }
        i += 1;
        if(i > 3)
            i = 0;
    }

    std::cout << count << std::endl;
}

int main(int argc, char * argv[]) {
    if(argc > 1) {
        std::ifstream infile(argv[1]);
        count_lines(infile);
    } else {
        count_lines(std::cin);
    }
}
