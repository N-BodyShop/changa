#pragma once

#include <algorithm>
#include <array>
#include <cstdio>
#include <string>
#include <vector>

/**
 * \class formatted_string
 *
 * \brief A simple formatted string to replace `s[n]printf` usage
 *
 * This is a simple handle class that provides storage for a formatted string
 * that would normally be allocated using a C-style array and populated with
 * `s[n]printf`. Storage is allocated either on the stack or the heap, depending
 * on the runtime size of the input as determined by calling
 * `std::snprintf(nullptr, 0, format, args);`. The underlying string can
 * be read using the class members `c_str` to retrieve a `char const*` or
 * `to_string` to retrieve a `std::string`.
 *
 * \warning Modifying the underlying data store is undefined behavior!
 *
 *
 * ===
 * Example usage:
 * @code
 *     // Old way
 *     char str[100];
 *     sprintf(str, "%s %d", "Hello", 3);
 *
 *     // New way
 *     auto str = make_formatted_string("%s %d", "Hello", 3);
 * @endcode
 */
template <typename... Args> class formatted_string {
public:
  /// Maximum length of the formatted string before heap allocation is used.
  static constexpr auto max_len = 32;

private:
  /*
   * The storage for the formatted string
   * A union is used to save space.
   */
  union arena {
    std::array<char, max_len> array;
    std::vector<char> vector;
    arena() {}
    ~arena() {}
  };
  arena storage;

  /// length of string (not including null-terminator)
  int len = 0;

  bool is_short() const noexcept { return len < max_len; }

public:
  /**
   * \param format The C-style printf format string
   * \param args A parameter pack of values to put into the string
   *
   * This uses the special form of `snprintf` that was introduced in C++11.
   * It determines the size of the resulting string that would be needed
   * to hold the formatted arguments. See
   * [cppref](https://en.cppreference.com/w/cpp/io/c/fprintf) for details.
   */
  explicit formatted_string(char const *format, Args... args) {
    static_assert(sizeof...(args) > 0, "Must specify at least one argument");

    len = std::snprintf(nullptr, 0, format, args...);
    char *const buffer = [&]() {
      if (is_short()) {
        new (&storage.array) decltype(storage.array)();
        return storage.array.data();
      } else {
        new (&storage.vector) decltype(storage.vector)(
            static_cast<size_t>(len) + 1);
        return storage.vector.data();
      }
    }();
    std::snprintf(buffer, static_cast<size_t>(len + 1), format, args...);
  }

  /// \brief Move constructor
  template <typename... Ts>
  formatted_string(formatted_string<Ts...> &&rhs) noexcept {
    this->len = rhs.len;
    if (is_short()) {
      new (&this->storage.array) decltype(this->storage.array)();
      std::copy(rhs.storage.array.begin(), rhs.storage.array.end(),
                this->storage.array.begin());
    } else {
      new (&this->storage.vector) decltype(storage.vector)(
          std::move(rhs.storage.vector));
    }
  }

  /// \brief Move assignment operator
  template <typename... Ts>
  formatted_string &operator=(formatted_string<Ts...> &&rhs) noexcept {
    new (this) formatted_string<Ts...>(std::move(rhs));
    return *this;
  }

  /**
   * Copy semantics make little sense, so they are explicitly deleted.
   */
  formatted_string(formatted_string const &) = delete;
  formatted_string &operator=(formatted_string const &) = delete;

  /**
   * Retrieve the underlying string as a `char const*`.
   * This is much like the `c_str` member of `std::string`.
   */
  char const *c_str() const noexcept {
    return (is_short()) ? storage.array.data() : storage.vector.data();
  }

  /**
   * Retrieve the underlying formatted string and place it in an
   * owning `std::string` object.
   */
  std::string to_string() const { return c_str(); }

  /**
   * The destructor only does work when the underlying string is long
   * enough to require heap storage (len > max_len).
   */
  ~formatted_string() {
    if (!is_short()) {
      storage.vector.~vector();
    }
  }
};

/**
 * C++ cannot deduce parameter types for constructors (this is fixed in
 * C++17; see below). As is common in the STL (e.g., std::make_pair), a
 * helper function is provided to overcome this limitation.
 *
 * ===
 * This allows the user to write
 * @code
 *     auto str = make_formatted_string("%s %d", "Hello", 3);
 * @endcode
 * instead of
 * @code
 *     auto str = formatted_string<(char&)[5], int>("%s %d", "Hello", 3);
 * @endcode
 */
template <typename... Args>
formatted_string<Args...> make_formatted_string(char const *format,
                                                Args... args) {
  return formatted_string<Args...>(format, args...);
}