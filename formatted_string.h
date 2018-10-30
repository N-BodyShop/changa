#pragma once

#include <string>
#include <array>
#include <vector>
#include <cstdio>
#include <algorithm>

template<typename... Args>
class formatted_string {
public:
    static constexpr auto max_len = 32;
private:
	union arena {
		std::array<char, max_len> array;
		std::vector<char> vector;
		arena() {}
		~arena() {}
	};
	arena storage;

	/// length of string (not including null-terminator)
    int len;

    bool is_short() const noexcept {
        return len < max_len;
    }
public:
    formatted_string(char const* format, Args... args) {
    	static_assert(sizeof...(args) > 0, "Must specify at least one argument");
		len = std::snprintf(nullptr, 0, format, args...);
		char *const buffer = [&]() {
			if(is_short()) {
				new (&storage.array) decltype(storage.array)();
				return storage.array.data();
			} else {
				new (&storage.vector) decltype(storage.vector)(static_cast<size_t>(len) + 1);
				return storage.vector.data();
			}
		}();
		std::snprintf(buffer, static_cast<size_t>(len + 1), format, args...);
    }

    template <typename... Ts>
    formatted_string(formatted_string<Ts...> &&rhs) {
    	this->len = rhs.len;
    	if (is_short()) {
    		new (&this->storage.array) decltype(this->storage.array)();
    		std::copy(rhs.storage.array.begin(), rhs.storage.array.end(), this->storage.array.begin());
    	} else {
    		new (&this->storage.vector) decltype(storage.vector)(std::move(rhs.storage.vector));
    	}
    }

    template <typename... Ts>
    formatted_string& operator=(formatted_string<Ts...> &&rhs) {
    	new (this) formatted_string<Ts...>(std::move(rhs));
    	return *this;
    }

    formatted_string(formatted_string const&) = delete;
    formatted_string& operator=(formatted_string const&) = delete;

    char const* c_str() const noexcept {
        return (is_short()) ? storage.array.data() : storage.vector.data();
    }
    std::string to_string() const {
        return c_str();
    }

    ~formatted_string() {
    	if(!is_short()) {
    		storage.vector.~vector();
    	}
    }
};

template<typename ... Args>
formatted_string<Args...> make_formatted_string(char const *format, Args ... args) {
    return formatted_string<Args...>(format, args...);
}
