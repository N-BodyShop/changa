#include "catch2.hpp"
#include "formatted_string.h"
#include <cstring>
#include <utility>
#include <string>

TEST_CASE( "string copy", "[formatted_string]" ) {
	char const* s = "Hello, world!\n";
	auto fs = make_formatted_string("%s", s);
	REQUIRE(std::strcmp(fs.c_str(), s) == 0);
}

TEST_CASE("move semantics", "[formatted_string]") {
	char const short_str[] = "Hello";
	char const long_str[] = "Lorem ipsum dolor sit amet, consectetur";

	static_assert(sizeof(short_str) < formatted_string<>::max_len);
	static_assert(sizeof(long_str) > formatted_string<>::max_len);

	auto const config = {std::make_pair("short",short_str),
						 std::make_pair("long", long_str)};

	for(auto const& c : config) {
		auto name = c.first;
		auto str = c.second;

		auto fs = make_formatted_string("%s", str);

		SECTION("move assign " + std::string(name)){
			auto fs2 = make_formatted_string("%d", 3);
			fs2 = std::move(fs);
			REQUIRE(std::strcmp(fs2.c_str(),str) == 0);
		}
		SECTION("move construct " + std::string(name)) {
			auto fs2 {std::move(fs)};
			REQUIRE(std::strcmp(fs2.c_str(),str) == 0);
		}
	}
}
