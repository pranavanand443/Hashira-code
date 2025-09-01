/**
 * Polynomial Solver - Shamir's Secret Sharing Implementation
 * 
 * This program solves polynomial coefficients using Shamir's Secret Sharing scheme.
 * It reads test cases in JSON format and uses Lagrange interpolation to find the 
 * constant term (secret) of the polynomial.
 * 
 * Author: GitHub Repository
 * Date: September 2025
 * Version: 2.0
 * 
 * Usage:
 *   ./polynomial_solver                          # Interactive mode with built-in test cases
 *   ./polynomial_solver < input.json            # Read JSON from stdin
 *   ./polynomial_solver input.json              # Read JSON from file
 *   ./polynomial_solver --test                  # Run comprehensive tests
 * 
 * Algorithm: Lagrange Interpolation
 * For a polynomial P(x) of degree m, given k = m + 1 points (x₁, y₁), ..., (xₖ, yₖ):
 * P(0) = Σᵢ₌₁ᵏ yᵢ * Πⱼ₌₁,ⱼ≠ᵢᵏ [(0 - xⱼ) / (xᵢ - xⱼ)]
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <climits>
#include <stdexcept>  // Added: For proper exception handling

using namespace std;

class PolynomialSolver {
private:
    struct Point {
        long long x;
        long double y;
        
        Point(long long x_val, long double y_val) : x(x_val), y(y_val) {}
    };

    /**
     * Convert a number from any base (2-16) to decimal
     * @param value: String representation of the number
     * @param base: Base of the number system (2-16)
     * @return: Decimal value as long double for precision
     * @throws invalid_argument: For invalid input
     */
    long double convertToDecimal(const string& value, int base) {
        if (value.empty() || base < 2 || base > 16) {
            throw invalid_argument("Invalid base (" + to_string(base) + ") or empty value");
        }
        
        long double result = 0.0L;
        long double power = 1.0L;
        
        // Process digits from right to left
        for (int i = value.length() - 1; i >= 0; i--) {
            char digit = tolower(value[i]);  // Fixed: Convert to lowercase for consistency
            int digitValue;
            
            if (digit >= '0' && digit <= '9') {
                digitValue = digit - '0';
            } else if (digit >= 'a' && digit <= 'f') {
                digitValue = digit - 'a' + 10;
            } else {
                throw invalid_argument("Invalid character '" + string(1, digit) + "' in number");
            }
            
            if (digitValue >= base) {
                throw invalid_argument("Digit " + to_string(digitValue) + " invalid for base " + to_string(base));
            }
            
            result += digitValue * power;
            power *= base;
        }
        
        return result;
    }

    /**
     * Lagrange interpolation to find polynomial value at x
     * @param points: Vector of points (x, y)
     * @param k: Number of points to use for interpolation
     * @param x: Point to evaluate polynomial at (default: 0 for secret)
     * @return: Polynomial value at x
     * @throws invalid_argument: For insufficient points or duplicate x values
     */
    long double lagrangeInterpolation(const vector<Point>& points, int k, long double x = 0.0L) {
        if (k <= 0 || k > (int)points.size()) {
            throw invalid_argument("Invalid k value: " + to_string(k));
        }
        
        // Check for duplicate x values
        for (int i = 0; i < k; i++) {
            for (int j = i + 1; j < k; j++) {
                if (points[i].x == points[j].x) {
                    throw invalid_argument("Duplicate x values found: " + to_string(points[i].x));
                }
            }
        }
        
        long double result = 0.0L;
        
        for (int i = 0; i < k; i++) {
            long double term = points[i].y;
            
            // Calculate Lagrange basis polynomial Li(x)
            for (int j = 0; j < k; j++) {
                if (i != j) {
                    long double denominator = (long double)(points[i].x - points[j].x);
                    if (abs(denominator) < 1e-15) {  // Fixed: Check for near-zero denominator
                        throw invalid_argument("Points too close together for stable interpolation");
                    }
                    term *= (x - points[j].x) / denominator;
                }
            }
            
            result += term;
        }
        
        return result;
    }

    /**
     * Simple JSON value extractor (no external dependencies)
     * @param json: JSON string
     * @param key: Key to extract
     * @return: Value as string
     */
    string extractValue(const string& json, const string& key) {
        string searchKey = "\"" + key + "\"";
        size_t keyPos = json.find(searchKey);
        if (keyPos == string::npos) return "";
        
        size_t colonPos = json.find(":", keyPos);
        if (colonPos == string::npos) return "";
        
        size_t valueStart = json.find("\"", colonPos);
        if (valueStart == string::npos) return "";
        valueStart++;
        
        size_t valueEnd = json.find("\"", valueStart);
        if (valueEnd == string::npos) return "";
        
        return json.substr(valueStart, valueEnd - valueStart);
    }
    
    /**
     * Extract integer value from JSON
     * @param json: JSON string
     * @param key: Key to extract
     * @return: Integer value or -1 if not found
     */
    int extractNumber(const string& json, const string& key) {
        string searchKey = "\"" + key + "\"";
        size_t keyPos = json.find(searchKey);
        if (keyPos == string::npos) return -1;
        
        size_t colonPos = json.find(":", keyPos);
        if (colonPos == string::npos) return -1;
        
        size_t numStart = colonPos + 1;
        while (numStart < json.length() && isspace(json[numStart])) {
            numStart++;
        }
        
        size_t numEnd = numStart;
        while (numEnd < json.length() && isdigit(json[numEnd])) {
            numEnd++;
        }
        
        if (numEnd > numStart) {
            try {
                return stoi(json.substr(numStart, numEnd - numStart));
            } catch (const exception&) {
                return -1;  // Fixed: Handle stoi exceptions
            }
        }
        
        return -1;
    }

public:
    /**
     * Solve polynomial from JSON input
     * @param jsonContent: JSON string containing the test case
     * @return: The secret (constant term) or LLONG_MIN on error
     */
    long long solveFromJSON(const string& jsonContent) {
        try {
            if (jsonContent.empty()) {
                cerr << "Error: Empty JSON content" << endl;
                return LLONG_MIN;
            }
            
            int n = extractNumber(jsonContent, "n");
            int k = extractNumber(jsonContent, "k");
            
            if (n <= 0 || k <= 0 || k > n) {  // Fixed: Added k > n check
                cerr << "Error: Invalid n=" << n << " or k=" << k << " (k must be ≤ n)" << endl;
                return LLONG_MIN;
            }
            
            cout << "Input: n=" << n << " roots, k=" << k << " minimum required" << endl;
            
            vector<Point> points;
            
            // Extract and convert all points
            for (int i = 1; i <= n; i++) {
                string pointKey = "\"" + to_string(i) + "\"";
                size_t pointStart = jsonContent.find(pointKey);
                if (pointStart == string::npos) continue;
                
                size_t braceStart = jsonContent.find("{", pointStart);
                size_t braceEnd = jsonContent.find("}", braceStart);
                if (braceStart == string::npos || braceEnd == string::npos) continue;
                
                string pointJson = jsonContent.substr(braceStart, braceEnd - braceStart + 1);
                
                string baseStr = extractValue(pointJson, "base");
                string valueStr = extractValue(pointJson, "value");
                
                if (!baseStr.empty() && !valueStr.empty()) {
                    try {
                        int base = stoi(baseStr);
                        long double decimalValue = convertToDecimal(valueStr, base);
                        points.push_back(Point(i, decimalValue));
                        
                        cout << "  Point " << i << ": \"" << valueStr << "\" (base " << base 
                             << ") = " << fixed << setprecision(0) << decimalValue << endl;
                    } catch (const exception& e) {
                        cerr << "  Warning: Skipping point " << i << " - " << e.what() << endl;
                        continue;
                    }
                }
            }
            
            if ((int)points.size() < k) {
                cerr << "Error: Not enough valid points (" << points.size() 
                     << " found, " << k << " required)" << endl;
                return LLONG_MIN;
            }
            
            // Use only the first k points for interpolation
            points.resize(k);
            
            // Use Lagrange interpolation to find the secret
            long double secret = lagrangeInterpolation(points, k, 0.0L);
            
            // Handle large numbers properly
            cout << "Secret (constant term): " << fixed << setprecision(0) << secret << endl;
            
            // Return as long long if it fits, otherwise return special value
            if (secret >= LLONG_MIN && secret <= LLONG_MAX) {
                return (long long)round(secret);
            } else {
                cout << "Note: Result exceeds long long range (returning 0)" << endl;
                return 0; // Fixed: Return 0 instead of -1 for overflow
            }
            
        } catch (const exception& e) {
            cerr << "Error processing JSON: " << e.what() << endl;
            return LLONG_MIN;
        }
    }

    /**
     * Run comprehensive tests
     */
    void runTests() {
        cout << "=== Running Comprehensive Tests ===" << endl;
        int passed = 0, total = 0;
        
        // Test 1: Base conversions
        cout << "\nTesting base conversions..." << endl;
        total++; if (abs(convertToDecimal("111", 2) - 7) < 0.01) { cout << "✓ Binary conversion"; passed++; } else cout << "✗ Binary conversion";
        total++; if (abs(convertToDecimal("213", 4) - 39) < 0.01) { cout << " ✓ Quaternary conversion"; passed++; } else cout << " ✗ Quaternary conversion";
        total++; if (abs(convertToDecimal("FF", 16) - 255) < 0.01) { cout << " ✓ Hex uppercase"; passed++; } else cout << " ✗ Hex uppercase";
        total++; if (abs(convertToDecimal("ff", 16) - 255) < 0.01) { cout << " ✓ Hex lowercase"; passed++; } else cout << " ✗ Hex lowercase";
        total++; if (abs(convertToDecimal("377", 8) - 255) < 0.01) { cout << " ✓ Octal conversion"; passed++; } else cout << " ✗ Octal conversion";
        cout << endl;
        
        // Test 2: Error handling
        cout << "\nTesting error handling..." << endl;
        total++;
        try {
            convertToDecimal("Z", 10);
            cout << "✗ Should catch invalid character";
        } catch (...) {
            cout << "✓ Catches invalid character";
            passed++;
        }
        
        total++;
        try {
            convertToDecimal("9", 8);
            cout << " ✗ Should catch invalid digit for base";
        } catch (...) {
            cout << " ✓ Catches invalid digit for base";
            passed++;
        }
        
        total++;
        try {
            convertToDecimal("", 10);
            cout << " ✗ Should catch empty string";
        } catch (...) {
            cout << " ✓ Catches empty string";
            passed++;
        }
        cout << endl;
        
        // Test 3: Known polynomial interpolation
        cout << "\nTesting polynomial interpolation..." << endl;
        vector<Point> testPoints = {Point(1, 1), Point(2, 4), Point(3, 9)}; // y = x^2
        long double result = lagrangeInterpolation(testPoints, 3, 0); // Should be 0
        total++;
        if (abs(result) < 0.01) {
            cout << "✓ Polynomial y=x² gives correct constant term (0)";
            passed++;
        } else {
            cout << "✗ Polynomial y=x² failed (got " << result << ")";
        }
        
        testPoints = {Point(1, 2), Point(2, 3), Point(3, 4)}; // y = x + 1
        result = lagrangeInterpolation(testPoints, 3, 0); // Should be 1
        total++;
        if (abs(result - 1.0) < 0.01) {
            cout << " ✓ Polynomial y=x+1 gives correct constant term (1)";
            passed++;
        } else {
            cout << " ✗ Polynomial y=x+1 failed (got " << result << ")";
        }
        
        // Test 4: Edge cases
        testPoints = {Point(0, 5), Point(1, 5), Point(2, 5)}; // y = 5 (constant)
        result = lagrangeInterpolation(testPoints, 3, 0); // Should be 5
        total++;
        if (abs(result - 5.0) < 0.01) {
            cout << " ✓ Constant polynomial y=5";
            passed++;
        } else {
            cout << " ✗ Constant polynomial failed (got " << result << ")";
        }
        cout << endl;
        
        // Test 5: Duplicate x values (should fail)
        cout << "\nTesting error conditions..." << endl;
        total++;
        try {
            testPoints = {Point(1, 1), Point(1, 2), Point(2, 3)};
            lagrangeInterpolation(testPoints, 3, 0);
            cout << "✗ Should catch duplicate x values";
        } catch (...) {
            cout << "✓ Catches duplicate x values";
            passed++;
        }
        cout << endl;
        
        cout << "Test Results: " << passed << "/" << total << " passed" << endl;
        if (passed == total) {
            cout << "🎉 All tests passed!" << endl;
        } else {
            cout << "⚠️  " << (total - passed) << " test(s) failed." << endl;
        }
    }

    /**
     * Get built-in test cases
     */
    vector<string> getTestCases() {
        return {
            // Test Case 1: Simple case with known answer
            R"({
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
            })",
            
            // Test Case 2: Complex case with large numbers
            R"({
                "keys": {
                    "n": 10,
                    "k": 7
                },
                "1": {
                    "base": "6",
                    "value": "13444211440455345511"
                },
                "2": {
                    "base": "15",
                    "value": "aed7015a346d635"
                },
                "3": {
                    "base": "15",
                    "value": "6aeeb69631c227c"
                },
                "4": {
                    "base": "16",
                    "value": "e1b5e05623d881f"
                },
                "5": {
                    "base": "8",
                    "value": "316034514573652620673"
                },
                "6": {
                    "base": "3",
                    "value": "2122212201122002221120200210011020220200"
                },
                "7": {
                    "base": "3",
                    "value": "20120221122211000100210021102001201112121"
                },
                "8": {
                    "base": "6",
                    "value": "20220554335330240002224253"
                },
                "9": {
                    "base": "12",
                    "value": "45153788322a1255483"
                },
                "10": {
                    "base": "7",
                    "value": "1101613130313526312514143"
                }
            })"
        };
    }
};

/**
 * Read entire file content
 * @param filename: Path to the file
 * @return: File content as string
 * @throws runtime_error: If file cannot be opened
 */
string readFile(const string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Cannot open file: " + filename);
    }
    
    ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

/**
 * Read from stdin
 * @return: stdin content as string
 */
string readStdin() {
    ostringstream ss;
    string line;
    while (getline(cin, line)) {
        ss << line << "\n";
    }
    return ss.str();
}

/**
 * Show usage information
 * @param programName: Name of the executable
 */
void showUsage(const char* programName) {
    cout << "Polynomial Solver - Shamir's Secret Sharing Implementation v2.0\n\n";
    cout << "Usage:\n";
    cout << "  " << programName << "                    # Interactive mode with built-in test cases\n";
    cout << "  " << programName << " --test            # Run comprehensive tests\n";
    cout << "  " << programName << " <file.json>       # Read JSON from file\n";
    cout << "  " << programName << " < input.json      # Read JSON from stdin\n";
    cout << "  " << programName << " --help            # Show this help\n\n";
    cout << "JSON Format:\n";
    cout << "{\n";
    cout << "  \"keys\": { \"n\": 4, \"k\": 3 },\n";
    cout << "  \"1\": { \"base\": \"10\", \"value\": \"4\" },\n";
    cout << "  \"2\": { \"base\": \"2\", \"value\": \"111\" },\n";
    cout << "  \"3\": { \"base\": \"10\", \"value\": \"12\" },\n";
    cout << "  \"6\": { \"base\": \"4\", \"value\": \"213\" }\n";
    cout << "}\n\n";
    cout << "Where:\n";
    cout << "  n = total number of roots provided\n";
    cout << "  k = minimum number of roots needed (polynomial degree + 1)\n";
    cout << "  base = number base (2-16)\n";
    cout << "  value = number in the specified base\n";
}

int main(int argc, char* argv[]) {
    try {
        PolynomialSolver solver;
        
        // Handle command line arguments
        if (argc > 1) {
            string arg = argv[1];
            
            if (arg == "--help" || arg == "-h") {
                showUsage(argv[0]);
                return 0;
            }
            
            if (arg == "--test") {
                solver.runTests();
                return 0;
            }
            
            if (arg == "--version" || arg == "-v") {
                cout << "Polynomial Solver v2.0" << endl;
                return 0;
            }
            
            // Try to read from file
            try {
                string content = readFile(arg);
                cout << "Reading from file: " << arg << endl;
                long long result = solver.solveFromJSON(content);
                if (result != LLONG_MIN) {
                    cout << "\nFinal Answer: " << result << endl;
                }
                return (result != LLONG_MIN) ? 0 : 1;
            } catch (const exception& e) {
                cerr << "Error reading file: " << e.what() << endl;
                return 1;
            }
        }
        
        // Check if stdin has data
        if (!cin.eof() && cin.peek() != EOF) {
            try {
                string content = readStdin();
                if (!content.empty()) {
                    cout << "Reading from stdin..." << endl;
                    long long result = solver.solveFromJSON(content);
                    if (result != LLONG_MIN) {
                        cout << "\nFinal Answer: " << result << endl;
                    }
                    return (result != LLONG_MIN) ? 0 : 1;
                }
            } catch (const exception& e) {
                cerr << "Error reading stdin: " << e.what() << endl;
                return 1;
            }
        }
        
        // Interactive mode with built-in test cases
        cout << "=== Polynomial Solver v2.0 - Interactive Mode ===" << endl;
        cout << "Running built-in test cases...\n" << endl;
        
        vector<string> testCases = solver.getTestCases();
        
        for (size_t i = 0; i < testCases.size(); i++) {
            cout << "--- Test Case " << (i + 1) << " ---" << endl;
            long long result = solver.solveFromJSON(testCases[i]);
            if (result != LLONG_MIN) {
                cout << "Final Answer: " << result << "\n" << endl;
            } else {
                cout << "Failed to solve this test case\n" << endl;
            }
        }
        
        cout << "Additional options:" << endl;
        cout << "  " << argv[0] << " --test     # Run comprehensive tests" << endl;
        cout << "  " << argv[0] << " --help     # Show detailed usage" << endl;
        cout << "  " << argv[0] << " file.json  # Process your own JSON file" << endl;
        
        return 0;
        
    } catch (const exception& e) {
        cerr << "Unexpected error: " << e.what() << endl;
        return 1;
    }
}
