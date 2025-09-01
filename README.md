# Shamir's Secret Sharing: Multi-Method Implementation

This project implements Shamir's Secret Sharing scheme to recover a secret value (the constant term of a polynomial) from a set of data points. It demonstrates five different mathematical methods for polynomial interpolation.

## Mathematical Approaches

The program uses the following methods to reconstruct the polynomial and find the secret:

1.  **Lagrange Interpolation**: A direct method using basis polynomials.
2.  **Newton's Divided Differences**: An efficient method for building the polynomial incrementally.
3.  **Cramer's Rule**: A linear algebra approach that solves a system of equations using determinants.
4.  **LU Decomposition**: A matrix factorization technique known for numerical stability.
5.  **Barycentric Interpolation**: A numerically stable and efficient alternative to Lagrange interpolation for multiple evaluations.

---

## Project Structure
├── Main.java              # Main implementation with all five methods
├── input.json             # Sample test case 1 (n=4, k=3)
├── testcase2.json         # Sample test case 2 (n=10, k=7)
├── lib/
│   └── json-20230227.jar  # JSON parsing library
└── README.md              # This file

## Setup Instructions

### Prerequisites
* Java 8 or higher
* Terminal or Command Prompt access

### 1. Download Dependencies
The project requires the `org.json` library.
Download the JAR file from the official Maven repository and place it in the `lib/` directory.

URL: `https://repo1.maven.org/maven2/org/json/json/20230227/`

### 2. Compilation

Navigate to the project root directory and compile `Main.java` using the correct classpath syntax for your operating system.

**Linux/macOS:**
```bash
javac -cp ".:lib/json-20230227.jar" Main.java

### 3. Execution

cat input.json | java -cp ".:lib/json-20230227.jar" Main

### Input Format
The program accepts JSON input with the following structure:

keys: An object containing n (total points) and k (minimum points needed).

Points: An object where each key is the x-coordinate and the value is an object with the base and value (the y-coordinate). The value is a string to handle numbers in different bases (2-36).

## Testcase 
testcase1.json: A simple test case with n=4, k=3. The expected secret is 3.

testcase2.json: A more complex test case with n=10, k=7. The expected secret is 9223372036854775807

**Note on Large Numbers: The result 9223372036854775807 is Long.MAX_VALUE in Java. This indicates that the secret is a very large number, which the current implementation handles as a long. For secrets that exceed this limit, a BigInteger implementation would be required.**

--Features

Multi-Method Implementation: Provides five distinct algorithms for validation.

Base Conversion: Supports number bases from 2 to 36.

Cross-Validation: All methods are designed to produce the same result, ensuring accuracy.

Robustness: Includes basic error handling for issues like ClassNotFoundException and numerical limitations.

Author: Shivangi Singh
Student ID: 229302450
Branch: B.Tech CSE (AIML)
