#pragma once

#include <cstring>
#include <string>
#include <vector>
#include <fstream>

namespace csv {
	class RowReader {
		std::ifstream ifstream;
		char columnSeparator = ',';
		char lineSeparator = '\n';
		char quote = '"';

	public:
		bool open(std::string filename) {
			ifstream = std::ifstream(filename.data(), std::ifstream::in);
			if (ifstream.good()) return true;
			return false;
		}
		
		bool done() {
			return !ifstream.good() || ifstream.eof();
		}
		
		std::vector<std::string> next() {
			std::vector<std::string> result;
			if (done()) return result;
			
			constexpr int N = 1024;
			char line[N];
			ifstream.getline(line, N, lineSeparator);
			// Split line into sections
			unsigned int length = strlen(line);
			unsigned int startIndex = 0, index = 0;
			while (index <= length) {
				if (index == length || line[index] == columnSeparator) {
					unsigned int endIndex = index;
					// Strip quotes
					if (line[startIndex] == quote) startIndex++;
					if (endIndex > 0 && line[endIndex - 1] == quote) endIndex--;

					std::string data(line + startIndex, endIndex - startIndex);
					result.push_back(data);

					index++;
					startIndex = index;
				} else {
					if (line[index] == quote) { // Skip ahead to ending quote
						index++;
						while (index < length && line[index] != quote) {
							index++;
						}
					}
					index++;
				}
			}
			return result;
		}
	};
}
