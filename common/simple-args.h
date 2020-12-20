#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <map>

#include "console-colours.h"

/** Expected use:

	SimpleArgs args(argc, argv);

	std::string foo = args.arg<std::string>("foo");
	std::string bar = args.arg<std::string>("bar", "a string for Bar", "default");

	// Exits if "foo" not supplied.  "bar" has a default value, so it's fine to omit
	if (args.help(std::cerr)) return 1;

**/
class SimpleArgs {
	int argc;
	char** argv;

	template<typename T>
	T valueFromString(const char *arg);
	
	std::string parsedCommand;
	struct Keywords {
		std::string keyword;
		std::string description;
		bool isHelp;
	};
	std::vector<Keywords> keywordOptions;
	std::vector<Keywords> argDetails;
	std::map<std::string, Keywords> flagOptions;
	void clearKeywords() {
		keywordOptions.resize(0);
		flagOptions.clear();
	}
	
	bool helpMode = false;
	bool hasError = false;
	std::string errorMessage;
	void setError(std::string message) {
		if (!hasError) {
			hasError = true;
			errorMessage = message;
		}
	}

	std::map<std::string, std::string> flagMap;
	void consumeFlags() {
		while (index < argc && std::strlen(argv[index]) > 0 && argv[index][0] == '-') {
			const char* arg = argv[index++];
			size_t length = strlen(arg);
			
			size_t keyStart = 1, keyEnd = keyStart + 1;
			size_t valueStart = keyEnd;
			// If it's "--long-arg" format
			if (length > 1 && arg[1] == '-') {
				keyStart++;
				while (keyEnd < length && arg[keyEnd] != '=') {
					keyEnd++;
				}
				valueStart = keyEnd;
				if (keyEnd < length) valueStart++;
			}

			std::string key = std::string(arg + keyStart, keyEnd - keyStart);
			std::string value = std::string(arg + valueStart);

			flagMap[key] = value;
	 	}
	 }
	
	int index = 1;
public:
	SimpleArgs(int argc, char* argv[]) : argc(argc), argv(argv) {
		parsedCommand = argv[0];
	}

	int help(std::ostream& out=std::cerr) {
		if (keywordOptions.size() > 0) {
			parsedCommand += std::string(" <command>");
		}
		out << "Usage:\n\t" <<  parsedCommand << "\n\n";
		if (keywordOptions.size() > 0) {
			out << "Commands:\n";
			for (unsigned int i = 0; i < keywordOptions.size(); i++) {
				out << "\t" << keywordOptions[i].keyword;
				if (keywordOptions[i].isHelp) out << " ...";
				if (keywordOptions[i].description.size()) out << "  -  " << keywordOptions[i].description;
				out << "\n";
			}
			out << "\n";
		}
		if (argDetails.size() > 0) {
			out << "Arguments:\n";
			for (auto iter = flagOptions.begin(); iter != flagOptions.end(); iter++) {
				Keywords &pair = iter->second;
				out << "\t" << (pair.keyword.length() > 1 ? "--" : "-") << pair.keyword;
				if (pair.description.size()) out << "  -  " << pair.description;
				out << "\n";
			}
			out << "\n";
		}
		if (flagOptions.size() > 0) {
			out << "Options:\n";
			for (auto iter = flagOptions.begin(); iter != flagOptions.end(); iter++) {
				Keywords &pair = iter->second;
				out << "\t" << (pair.keyword.length() > 1 ? "--" : "-") << pair.keyword;
				if (pair.description.size()) out << "  -  " << pair.description;
				out << "\n";
			}
			out << "\n";
		}
		return hasError ? -1 : 0;
	}

	bool error(std::ostream& out=std::cerr) {
		if (!hasError && !helpMode) return false;
		help(out);
		if (!helpMode) {
			out << Console::Red << errorMessage << Console::Reset << "\n";
		}
		return true;
	}
	bool error(std::string forcedError, std::ostream& out=std::cerr) {
		if (!hasError) {
			hasError = true;
			errorMessage = forcedError;
		}
		return error(out);
	}
	
	template<typename T>
	T arg(std::string name, std::string longName, T defaultValue) {
		clearKeywords();
		consumeFlags();
		parsedCommand += std::string(" [?") + name + "]";
		argDetails.push_back(Keywords{name, longName, false});

		if (index >= argc) return defaultValue;
		return valueFromString<T>(argv[index++]);
	}

	template<typename T>
	T arg(std::string name, std::string longName="") {
		clearKeywords();
		consumeFlags();
		parsedCommand += std::string(" [") + name + "]";
		argDetails.push_back(Keywords{name, longName, false});

		if (index >= argc) {
			if (longName.length() > 0) {
				setError("Missing " + longName + " <" + name + ">");
			} else {
				setError("Missing argument <" + name + ">");
			}
			return T();
		}

		return valueFromString<T>(argv[index++]);
	}

	bool command(std::string keyword, std::string description="", bool isHelp=false) {
		consumeFlags();
		if (index < argc && !keyword.compare(argv[index])) {
			clearKeywords();
			index++;
			if (!isHelp) parsedCommand += std::string(" ") + keyword;
			return true;
		}
		keywordOptions.push_back(Keywords{keyword, description, isHelp});
		return false;
	}
	bool helpCommand(std::string keyword) {
		helpMode = command(keyword, "", true);
		if (helpMode) {
			keywordOptions.insert(keywordOptions.begin(), Keywords{keyword, "", true});
		}
		return helpMode;
	}

	template<typename T=std::string>
	T flag(std::string key, std::string description, T defaultValue) {
		consumeFlags();
		if (!hasFlag(key, description)) return defaultValue;

		auto iterator = flagMap.find(key);
		return valueFromString<T>(iterator->second.c_str());
	}
	template<typename T=std::string>
	T flag(std::string key, T defaultValue) {
		consumeFlags();
		if (!hasFlag(key, "")) return defaultValue;

		auto iterator = flagMap.find(key);
		return valueFromString<T>(iterator->second.c_str());
	}
	template<typename T=std::string>
	T flag(std::string key) {
		return flag<T>(key, T());
	}
	bool hasFlag(std::string key, std::string description="") {
		consumeFlags();
		auto iterator = flagMap.find(key);
		if (description.length() > 0 || iterator == flagMap.end()) {
			flagOptions[key] = Keywords{key, description, false};
		}

		iterator = flagMap.find(key);
		return iterator != flagMap.end();
	}
	bool helpFlag(std::string key, std::string description="") {
		consumeFlags();
		flagOptions[key] = Keywords{key, description, true};
		auto iterator = flagMap.find(key);
		helpMode = (iterator != flagMap.end());
		return helpMode;
	}
};

template<>
std::string SimpleArgs::valueFromString(const char *arg) {
	return arg;
}
template<>
const char * SimpleArgs::valueFromString(const char *arg) {
	return arg;
}
template<>
int SimpleArgs::valueFromString(const char *arg) {
	return std::stoi(arg);
}
template<>
long SimpleArgs::valueFromString(const char *arg) {
	return std::stol(arg);
}
template<>
float SimpleArgs::valueFromString(const char *arg) {
	return std::stof(arg);
}
template<>
double SimpleArgs::valueFromString(const char *arg) {
	return std::stod(arg);
}