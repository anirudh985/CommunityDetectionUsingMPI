#include "ReformatFile.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>

using namespace std;

/*
*	Maximum word length for wiki-Talk.
*   Taking it as 7, as the number of edges is 7 digit number
*/

int lineSize = 0;
string appendZero(string s, int currLen);
string appendZeroForDouble(string s, int currLen);
int getLineSize();

//int main(int argc, char* argv[]){
//	runFormattor();
//	return 0;
//}

void runFormattor() {
	ifstream infile;
//	infile.open("Text.txt");
	infile.open("soc-livejournal1.mtx");

	ofstream outfile;
//	outfile.open("out.txt", std::ofstream::out | std::ofstream::trunc);
	outfile.open("output", std::ofstream::out | std::ofstream::trunc);

	if (!outfile.is_open()) {
		cout << "Unable to open output file \n";
	}

	if (infile.is_open()) {		
		string line;
		while (getline(infile, line)) {
			istringstream iss(line);
			string s0, s1, s2, s3;
			
//			if (!(iss >> s0 >> s1 >> s2 >> s3)) {
//				cout << "Some weird error\n";
//			}

			if (!(iss >> s0 >> s1)) {
				cout << "Some weird error\n";
			}
			
			// Ignore s0
			s0 = appendZero(s0, s0.length());
			s1 = appendZero(s1, s1.length());
//			s3 = appendZeroForDouble(s3, s3.length());

//			outfile << s1 << "	" << s2 << "	"<<s3<<"\n";

			outfile << s0 << "	" << s1 <<"\n";
			outfile << s1 << "	" << s0 <<"\n";

		}

		// Line size will be equal to the length, because sizeOf(char) = 1 minus tab and \n
		lineSize = 4 * MAX_WORD_SIZE + 4;
		cout << "Line Size " << lineSize;
		cout << "Exiting Loop \n";
	} else {
		cout << "Unable to open file\n";
	}

	outfile.close();
	infile.close();
}

string appendZero(string s, int currLen) {
	for (int i = 0; i < MAX_WORD_SIZE - currLen; i++) {
		s = "0" + s;
	}
	return s;
}

string appendZeroForDouble(string s, int currLen) {
	for (int i = 0; i < 2 * MAX_WORD_SIZE - currLen; i++) {
		s = "0" + s;
	}
	return s;
}

int getLineSize() {
//	if (lineSize <= 0) {
//		return 4 * MAX_WORD_SIZE + 3;
//	}
//	return lineSize;
	if (lineSize <= 0) {
		return 3 * MAX_WORD_SIZE + 3;
	}
	return lineSize;
}
