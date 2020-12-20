//extern C {
#include <sqlite3.h>
//}

#include <memory>
#include <string>

#include "console-colours.h"

/* Expected use:

	SQLite3Versioned db("database.sqlite");
	// Backwards-compatible structure - append to this as changes are made
	db.bumpVersion(1, "CREATE TABLE hello_world (foo INTEGER, bar TEXT)");
	
	// Template-based query
	auto find = db.statement("SELECT bar FROM hello_world WHERE foo=@foo");
	// Specify values by index, or by string
	find->bind(1, 123);
	find->bind<double>("foo", 456); // can explicitly set type
	
	// Step through result
	if (find->step()) {
		std::string bar = find->get<std::string>(0);
		...
	}
	
	// Transactions are mostly used for efficiency.  If you have multiple threads, re-consider using SQLite.
	{
		auto transaction = db.transaction();
		
		auto query = db.statement(...);
	}
*/

class SQLite3 {
	std::string filename;
	sqlite3* m_db = nullptr;

	sqlite3* db() { // Lazy-loading
		if (!m_db) {
#ifdef DEBUG_SQL
	std::cout << Console::Blue << "lazy-loading SQL file: " << filename << Console::Reset << "\n";
#endif
			int resultCode = sqlite3_open(filename.c_str(), &m_db);
			if (resultCode) {
				error = 0;
				errorMessage = std::string(sqlite3_errmsg(m_db));
				sqlite3_close(m_db);
				m_db = nullptr;

				std::cerr << Console::Bright << Console::Red << "SQL error: " << error  << " - " << errorMessage << Console::Reset << "\n";
			}
		}
		return m_db;
	}
public:
	int error = 0;
	std::string errorMessage;

	SQLite3(std::string filename) : filename(filename) {}
	~SQLite3() {
		if (m_db) sqlite3_close(m_db);
	}
	
	class Statement {
		SQLite3* db;
		sqlite3_stmt* statement = 0;
		const char * sqlRemainder = 0;

		void getError(int code) {
			error = db->error = code;
			errorMessage = db->errorMessage = std::string(sqlite3_errmsg(db->db()));
		}
		void checkForError(int resultCode) {
			if (resultCode != SQLITE_OK) {
				getError(resultCode);
				//statement = nullptr;

				std::cerr << Console::Bright << Console::Red << "SQL error: " << error  << " - " << errorMessage << Console::Reset << "\n";
			}
		}
		void customError(std::string message) {
			error = -1;
			errorMessage = message;
			std::cerr << Console::Bright << Console::Red << "Error: " << message << Console::Reset << "\n";
		}
	public:
		int error = 0;
		std::string errorMessage;

		Statement(SQLite3* db, std::string sqlString) : db(db) {
#ifdef DEBUG_SQL
			std::cout << Console::Dim << sqlString << Console::Reset << "\n";
#endif
			int resultCode = sqlite3_prepare_v2(db->db(), sqlString.c_str(), sqlString.length(), &statement, &sqlRemainder);
			if (resultCode != SQLITE_OK) {
				getError(resultCode);
				statement = nullptr;

				std::cerr << Console::Bright << Console::Red << "SQL error: " << error  << " - " << errorMessage << Console::Reset << "\n";
			}
		}
		~Statement() {
			if (statement) sqlite3_finalize(statement);
		}
		
		template<typename T>
		void bind(const char* name, T value) {
			int index = sqlite3_bind_parameter_index(statement, name);
			this->bind(index, value);
		}
		
		void bind(int index, std::string value) {
			checkForError(sqlite3_bind_text(statement, index, value.data(), value.size(), SQLITE_TRANSIENT));
		}
		void bind(int index, char const *value) {
			checkForError(sqlite3_bind_text(statement, index, value, -1, SQLITE_TRANSIENT));
		}
		void bind(int index, double value) {
			checkForError(sqlite3_bind_double(statement, index, value));
		}
		void bind(int index, int value) {
			checkForError(sqlite3_bind_int(statement, index, value));
		}
		void bind(int index, long value) {
			checkForError(sqlite3_bind_int64(statement, index, value));
		}
		void bind(int index) {
			checkForError(sqlite3_bind_null(statement, index));
		}
		void reset() {
			checkForError(sqlite3_reset(statement));
		}

		bool step() {
			if (!statement) return false;
			int queryStatus = sqlite3_step(statement);
			if (queryStatus != SQLITE_ROW && queryStatus != SQLITE_DONE) {
				getError(queryStatus);
			}
			return queryStatus == SQLITE_ROW;
		}

		int columns() {
			return statement ? sqlite3_column_count(statement) : 0;
		}

		int columnType(int i) {
			return statement ? sqlite3_column_type(statement, i) : 0;
		}

		std::string columnName(int i) {
			if (!statement) return "";
			const char* cstr = sqlite3_column_name(statement, i);
			return std::string(reinterpret_cast<const char*>(cstr));
		}

		template <typename T=std::string>
		T get(int i);

		template <typename T=int>
		void bind(int i, T value);
	};

	std::shared_ptr<Statement> statement(std::string sqlString) {
		return std::make_shared<Statement>(this, sqlString);
	}
	std::shared_ptr<Statement> query(std::string sqlString) {
		return this->statement(sqlString);
	}

	sqlite3_int64 insertId() {
		return sqlite3_last_insert_rowid(db());
	}

	class Transaction {
		SQLite3* db;
		Statement open;
		Statement close;
	public:
		Transaction(SQLite3* db) : db(db), open(db, "BEGIN TRANSACTION"), close(db, "COMMIT TRANSACTION") {
			open.step();
		}
		~Transaction() {
			close.step();
			db->currentTransaction = nullptr;
		}
	};

	std::shared_ptr<Transaction> transaction() {
		if (currentTransaction) return nullptr;
		return currentTransaction = std::make_shared<Transaction>(this);
	}
private:
	std::shared_ptr<Transaction> currentTransaction;
};

template<>
int SQLite3::Statement::get<int>(int i) {
	return statement ? sqlite3_column_int(statement, i) : 0;
}
template<>
long SQLite3::Statement::get<long>(int i) {
	return statement ? sqlite3_column_int64(statement, i) : 0;
}
template<>
double SQLite3::Statement::get<double>(int i) {
	return statement ? sqlite3_column_double(statement, i) : 0;
}
template<>
std::string SQLite3::Statement::get<std::string>(int i) {
	if (!statement) return "";
	const unsigned char* cstr = sqlite3_column_text(statement, i);
	if (!cstr) return "";
	return std::string(reinterpret_cast<const char*>(cstr));
}

template<>
void SQLite3::Statement::bind<int>(int i, int value) {
	if (statement) sqlite3_bind_int(statement, i, value);
}

/********/

class SQLite3Versioned : public SQLite3 {
	int userVersion = 0;
	bool versionSuccess = true;
public:
	SQLite3Versioned(std::string filename) : SQLite3(filename) {
		auto query = statement("PRAGMA user_version");
		query->step();
		userVersion = query->get<int>(0);
	}

	/* Use to create sequence of changes which upgrade old DBs, like:
		db.bumpVersion(1, "CREATE TABLE hello_world (foo INTEGER, bar TEXT)");
		db.bumpVersion(2, "ALTER TABLE hello_world ADD beep FLOAT");
	*/
	bool bumpVersion(int version, std::string sql) {
		if (!versionSuccess) return false; // Fail once, don't continue
		if (userVersion < version) {
			userVersion = version;
			auto query = statement(sql);
			if (query->error) return versionSuccess = false;
			query->step();
			if (query->error) return versionSuccess = false;

			query = statement(std::string("PRAGMA user_version=") + std::to_string(version));
			if (query->error) return versionSuccess = false;
			query->step();
			return true;
		}
		return versionSuccess;
	}
};
