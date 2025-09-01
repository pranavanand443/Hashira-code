# Polynomial Solver - Shamir's Secret Sharing

A C++ implementation of Shamir's Secret Sharing scheme using Lagrange interpolation to find polynomial coefficients and recover secrets.

## 🎯 Overview

This program solves polynomial interpolation problems commonly found in cryptographic applications, specifically Shamir's Secret Sharing. Given a set of points (x, y) where y-values are encoded in various number bases, it reconstructs the original polynomial and finds the secret (constant term).

## ✨ Features

- **Multiple Input Methods**: JSON files, stdin, or built-in test cases
- **Base Conversion**: Supports bases 2-16 (binary, octal, decimal, hexadecimal, etc.)
- **High Precision**: Uses `long double` for accurate calculations with large numbers
- **Comprehensive Testing**: Built-in test suite with edge case validation
- **Self-Contained**: No external dependencies required
- **Cross-Platform**: Works on Windows, Linux, and macOS

## 🚀 Quick Start

### Compilation
```bash
g++ -std=c++17 -O2 -o polynomial_solver polynomial_solver.cpp
```

### Usage
```bash
# Interactive mode with built-in test cases
./polynomial_solver

# Run comprehensive tests
./polynomial_solver --test

# Process JSON from file
./polynomial_solver input.json

# Process JSON from stdin
./polynomial_solver < input.json
cat input.json | ./polynomial_solver

# Show help
./polynomial_solver --help
```

## 📋 Input Format

JSON format with the following structure:
```json
{
    "keys": {
        "n": 4,
        "k": 3
    },
    "1": {
        "base": "10",
        "value": "4"
    },
    "2": {
        "base": "2", 
        "value": "111"
    },
    "3": {
        "base": "10",
        "value": "12"
    },
    "6": {
        "base": "4",
        "value": "213"
    }
}
```

Where:
- `n`: Total number of roots provided
- `k`: Minimum number of roots required (k = polynomial_degree + 1)
- Each numbered key represents a point with its base and value

## 🧮 Algorithm

The program uses **Lagrange Interpolation** to reconstruct polynomials:

```
P(x) = Σᵢ₌₁ᵏ yᵢ × Πⱼ₌₁,ⱼ≠ᵢᵏ [(x - xⱼ) / (xᵢ - xⱼ)]
```

To find the secret (constant term), it evaluates P(0):

```
Secret = P(0) = Σᵢ₌₁ᵏ yᵢ × Πⱼ₌₁,ⱼ≠ᵢᵏ [(0 - xⱼ) / (xᵢ - xⱼ)]
```

## 📊 Test Cases

### Test Case 1
- **Input**: 4 points, need 3 minimum
- **Points**: (1,4), (2,7), (3,12), (6,39)
- **Result**: Secret = **3**

### Test Case 2  
- **Input**: 10 points, need 7 minimum
- **Complex**: Various bases (3, 6, 7, 8, 12, 15, 16)
- **Result**: Secret = **66983859479598506132**

## 🔍 Testing

The program includes comprehensive testing:
- Base conversion validation
- Error handling verification  
- Polynomial interpolation accuracy
- Edge case handling

```bash
./polynomial_solver --test
```

## 🏗️ Implementation Details

### Key Components
1. **Base Conversion**: Converts numbers from any base (2-16) to decimal
2. **JSON Parsing**: Lightweight parser without external dependencies
3. **Lagrange Interpolation**: Numerically stable algorithm for polynomial reconstruction
4. **Error Handling**: Robust validation and error reporting

### Why Lagrange Interpolation?
- More numerically stable than matrix methods for this specific problem
- Direct computation of the secret without needing full polynomial coefficients
- Fewer floating-point operations reduce accumulation errors
- Optimal for Shamir's Secret Sharing applications

## 🎓 Mathematical Background

**Shamir's Secret Sharing** is a cryptographic algorithm that divides a secret into multiple parts (shares). The secret can only be reconstructed when a sufficient number of shares are combined.

- **Security**: Based on polynomial interpolation over finite fields
- **Threshold**: Requires exactly `k` shares to reconstruct the secret
- **Applications**: Secure key storage, distributed systems, cryptographic protocols

## 🤝 Contributing

Feel free to contribute by:
- Adding support for more number bases
- Implementing finite field arithmetic
- Adding more test cases
- Improving documentation

## 📄 License

This project is open source and available under the MIT License.

## 🔗 References

- [Shamir's Secret Sharing](https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing)
- [Lagrange Interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial)
- [Polynomial Interpolation](https://mathworld.wolfram.com/PolynomialInterpolation.html)
