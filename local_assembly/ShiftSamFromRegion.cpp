#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string.h>

using namespace std;

int main(int argc, const char * argv[]) {
	string line;
	int index = 0;
	if (argc > 1 and (strcmp(argv[1], "-h") == 0  or  strcmp(argv[1], "--help") == 0)) {
		cout << endl 
				 << "Usage: cat file.sam | shiftSamPos > output.sam" << endl
				 << "  The reference name must be in the format chrom:start-end " << endl
				 << "  which is the default output from samtools faidx." << endl;
		exit(0);
	}
	bool preserveTitle = false;
	if (argc > 1 and strcmp(argv[1], "--preserveTitle") == 0) {
		preserveTitle = true;
	}

	while (getline(cin, line)) {
		if (line == "") {
			break;
		}
		if (line[0] == '@') {
			continue;
		}

		stringstream lineStrm(line);
		string val;
		vector<string> vals;
		while ( (lineStrm >> val) ) {
			if (val == "") {
				break;
			}
			vals.push_back(val);
		}

		stringstream name;
		name << vals[2] << "/" << index;
		if (preserveTitle == false) {
			vals[0] = name.str();
		}
		
		if (vals[2] != "*") {
			int i;
			i = 0;
			while (i < vals[2].size() && (vals[2][i] != ':' and vals[2][i] != '.')) { i++;}
			string chrom = vals[2].substr(0,i);
			i++;
			int posStart = i;
			while (i < vals[2].size() && vals[2][i] != '-') { i++; }
			string startString = vals[2].substr(posStart, i-posStart);
			i++;
			string endString  = vals[2].substr(i, vals[2].size() - i);
			int start = atoi(startString.c_str());
			int offset = atoi(vals[3].c_str());
			stringstream posStrm;
			posStrm << start + offset - 1;
			vals[2] = chrom;
			vals[3] = posStrm.str();
		}
		int i;
		for (i =0 ; i < vals.size(); i++) {
			cout << vals[i];
			if (i < vals.size()-1) {
				cout << "\t";
			}
		}
		cout << endl;
		index +=1;
	}
}
