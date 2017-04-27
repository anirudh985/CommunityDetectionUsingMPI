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

void runFormattor() {
	ifstream infile;
	infile.open("Text.txt");

	ofstream outfile;
	outfile.open("out.txt", std::ofstream::out | std::ofstream::trunc);

	if (!outfile.is_open()) {
		cout << "Unable to open output file \n";
	}

	if (infile.is_open()) {		
		string line;
		while (getline(infile, line)) {
			istringstream iss(line);
			string s0, s1, s2, s3;
			
			if (!(iss >> s0 >> s1 >> s2 >> s3)) {
				cout << "Some weird error\n";
			}
			
			// Ignore s0
			s1 = appendZero(s1, s1.length());
			s2 = appendZero(s2, s2.length());
			s3 = appendZeroForDouble(s3, s3.length());

			outfile << s1 << "	" << s2 << "	"<<s3<<"\n";

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
	if (lineSize <= 0) {
		return 4 * MAX_WORD_SIZE + 3;
	}
	return lineSize;
}