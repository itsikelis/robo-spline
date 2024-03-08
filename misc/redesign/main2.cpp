#include <iostream>

template <size_t _Dim>
class Foo {
private:
    std::array<double, Dim> _arr;

public:
    Foo() {}
};

int main()
{
    std::cout << "Hello from main2!" << std::endl;
    return 0;

    Foo<2> foo = Foo<2>();
}