/*
 * preprocess.cpp
 *
 *  Created on: May 4, 2017
 *      Author: osu8229
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>

#define MAX_WORD_SIZE 7

using namespace std;

string appendZero(string s, int currLen);
int getRenumberedVertex(unsigned long originalVertex, vector<unsigned long> &v);
void processOutput(vector<unsigned long> &v);
unsigned long binarySearch(unsigned long* arr, unsigned long value, unsigned long left, unsigned long right);
void insertIntoMap(unsigned long key, unsigned long value, map<unsigned long, vector<pair<unsigned long, long> > > &m);
vector<pair<unsigned long, long> >::iterator find1(vector<pair<unsigned long, long> >::iterator vBegin, vector<pair<unsigned long, long> >::iterator vEnd, unsigned long value);
int main(int argc, char* argv[]){

	set<unsigned long> vertices;
	vector<unsigned long> renumberedIndexVertices;
	unsigned long val;
	ifstream infile;
	infile.open("soc-livejournal1.mtx");
//	infile.open("dummy");

	if (infile.is_open()) {
		string line;
		while (getline(infile, line)) {
			istringstream iss(line);
			string s0, s1;

			if (!(iss >> s0 >> s1)) {
				cout << "Some weird error\n";
			}

			stringstream(s0) >> val;
			vertices.insert(val);

			stringstream(s1) >> val;
			vertices.insert(val);

		}

		cout << "Exiting Loop \n";
	} else {
		cout << "Unable to open file\n";
	}

	infile.close();

	vector<unsigned long> v(vertices.begin(), vertices.end());

	processOutput(v);

	return 0;
}


int getRenumberedVertex(unsigned long originalVertex, vector<unsigned long> &v){
//	return lower_bound(v.begin(), v.end(), originalVertex) - v.begin();
	return binarySearch(&v[0], originalVertex, 0, v.size());
}

string appendZero(string s, int currLen) {
	for (int i = 0; i < MAX_WORD_SIZE - currLen; i++) {
		s = "0" + s;
	}
	return s;
}

void processOutput(vector<unsigned long>& v){
	ifstream infile;
//	infile.open("soc-livejournal1.mtx");
//	infile.open("dummy");
	infile.open("co-papers-dblp.mtx");

	ofstream outfile;
	outfile.open("inputToLouvain_copapers", std::ofstream::out | std::ofstream::trunc);

	map<unsigned long, vector<pair<unsigned long, long> > > m;
	unsigned long key;
	unsigned long value;

	unsigned long dummyVal;

	if (infile.is_open()) {
		string line;
		while (getline(infile, line)) {
			istringstream iss(line);
			string s0, s1;

			if (!(iss >> s0 >> s1)) {
				cout << "Some weird error\n";
			}

			stringstream(s0) >> dummyVal;
			key = binarySearch(&v[0], dummyVal, 0, v.size());

			stringstream(s1) >> dummyVal;
			value = binarySearch(&v[0], dummyVal, 0, v.size());

			insertIntoMap(key, value, m);
			insertIntoMap(value, key, m);
		}
		cout<< "Exiting loop" << endl;
	} else {
		cout << "Unable to open file\n";
	}

	cout << "Beginning to write to outfile" << endl;
	cout << "size of map is "<< m.size() << endl;

	unsigned long i = 0;
	map<unsigned long, vector<pair<unsigned long, long> > >::iterator it;

	for(it = m.begin(); it != m.end(); it++){
		string s1;
		stringstream x;
		x << it->first;
		x >> s1;

		s1 = appendZero(s1, s1.length());

		vector<pair<unsigned long ,long> >::iterator vit;

		for(vit = it->second.begin(); vit != it->second.end(); vit++){
			stringstream y, z;
			string s2, s3;
			y << vit->first;
			y >> s2;
			s2 = appendZero(s2, s2.length());

			z << vit->second;
			z >> s3;
			s3 = appendZero(s3, s3.length());
			outfile << s1 << " " << s2 << " " << s3 <<"\n";
		}

		i++;
		if(i == m.size() / 2)
			cout << "Completed half" << endl;

	}


}

void insertIntoMap(unsigned long key, unsigned long value, map<unsigned long, vector<pair<unsigned long, long> > > &m){

	map<unsigned long, vector<pair<unsigned long, long> > >::iterator it;

	it = m.find(key);

	if(it != m.end()){	// found key in the map
		vector<pair<unsigned long, long> >::iterator vit;
		vit = find1(it->second.begin(), it->second.end(), value);
		if(vit == it->second.end()){ // didn't find destination vertex
			it->second.push_back(make_pair(value, 1));
		}
		else{
			(*vit).second++;
		}
	}
	else{				// didn't find the key
		m[key].push_back(make_pair(value, 1));
	}

	return;
}

unsigned long binarySearch(unsigned long* arr, unsigned long value, unsigned long left, unsigned long right){
      while (left <= right) {
            int middle = (left + right) / 2;
            if (arr[middle] == value)
                  return middle;
            else if (arr[middle] > value)
                  right = middle - 1;
            else
                  left = middle + 1;
      }
      return -1;
}

vector<pair<unsigned long, long> >::iterator find1(vector<pair<unsigned long, long> >::iterator vBegin, vector<pair<unsigned long, long> >::iterator vEnd, unsigned long value){

	vector<pair<unsigned long, long> >::iterator vit;

	for(vit = vBegin; vit != vEnd; vit++){
		if((*vit).first == value) break;
	}

	return vit;
}
