# Coding conventions

## general
* user defined types and type aliases use CamelCase,
  e.g., `DerivativeMethod` or `Interpolate`
* functions and variables use `lower_case` names, e.g.,
  `long_function_name`, `variable_name`, `derivative`, `sample_2d`
* 4 spaces are used instead of tabs
* whenever feasible, lines are at most 80 characters long
  
## c++ specific
* try to follow the
  [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
* the curly braces in function definitions are placed as in
  ```
  int my_function(int argument)
  {
      return argument * 2;
  }
  ```
  while everywhere else they are placed as in
  ```
  class ClassName {
    // stuff
  };
  ```
  or
  ```a
  if condition {
    // stuff
  }
  ```

## python specific
* [PEP8](https://www.python.org/dev/peps/pep-0008/) should be followed wherever
  feasible
* global constants are capitalized, e.g., `MY_CONSTANT`
* string are delimited using `'` instead of `"`