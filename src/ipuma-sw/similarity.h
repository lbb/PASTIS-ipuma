#ifndef SIMILARITY_H
#define SIMILARITY_H
#include<cstdint>
#include "matrix.h"
// Implementations for different similarity scores between single chars
// Will be extended in the future with specific scoresets for DNA/RNA or Protein.

namespace swatlib {

/**
 * Similarity score +1 if similar, -1 if different.
 */
int simpleSimilarity(char a, char b);

// NA Substituion matrix
const std::vector<int8_t> NA_SUBSTITUTION = {
//N  A  T  C  G  U  R  Y  K  M  S  W  B  D  H  V
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, // N
	1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, // A
	1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, // T
	1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, // C
	1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, // G
	1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // U
	1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, // R
	1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, // Y
	1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, // K
	1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, // M
	1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, // S
	1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, // W
	1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, // B
	1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, // D
	1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, // H
	1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, // V
};

// BLOSUM matrices taken from https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/
const std::vector<int8_t> BLOSUM50 = {
// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
   5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0,-2,-2,-1,-1,-5, // A
  -2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3,-1,-3, 0,-1,-5, // R
  -1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3, 5,-4, 0,-1,-5, // N
  -2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4, 6,-4, 1,-1,-5, // D
  -1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-3,-2,-3,-1,-5, // C
  -1, 1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3, 0,-3, 4,-1,-5, // Q
  -1, 0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3, 1,-3, 5,-1,-5, // E
   0,-3, 0,-1,-3,-2,-3, 8,-2,-4,-4,-2,-3,-4,-2, 0,-2,-3,-3,-4,-1,-4,-2,-1,-5, // G
  -2, 0, 1,-1,-3, 1, 0,-2,10,-4,-3, 0,-1,-1,-2,-1,-2,-3, 2,-4, 0,-3, 0,-1,-5, // H
  -1,-4,-3,-4,-2,-3,-4,-4,-4, 5, 2,-3, 2, 0,-3,-3,-1,-3,-1, 4,-4, 4,-3,-1,-5, // I
  -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5,-3, 3, 1,-4,-3,-1,-2,-1, 1,-4, 4,-3,-1,-5, // L
  -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6,-2,-4,-1, 0,-1,-3,-2,-3, 0,-3, 1,-1,-5, // K
  -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7, 0,-3,-2,-1,-1, 0, 1,-3, 2,-1,-1,-5, // M
  -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8,-4,-3,-2, 1, 4,-1,-4, 1,-4,-1,-5, // F
  -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3,-2,-3,-1,-1,-5, // P
   1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5, 2,-4,-2,-2, 0,-3, 0,-1,-5, // S
   0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5,-3,-2, 0, 0,-1,-1,-1,-5, // T
  -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15, 2,-3,-5,-2,-2,-1,-5, // W
  -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8,-1,-3,-1,-2,-1,-5, // Y
   0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5,-3, 2,-3,-1,-5, // V
  -2,-1, 5, 6,-3, 0, 1,-1, 0,-4,-4, 0,-3,-4,-2, 0, 0,-5,-3,-3, 6,-4, 1,-1,-5, // B
  -2,-3,-4,-4,-2,-3,-3,-4,-3, 4, 4,-3, 2, 1,-3,-3,-1,-2,-1, 2,-4, 4,-3,-1,-5, // J
  -1, 0, 0, 1,-3, 4, 5,-2, 0,-3,-3, 1,-1,-4,-1, 0,-1,-2,-2,-3, 1,-3, 5,-1,-5, // Z
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-5, // X
  -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5, 1, // *
};

const std::vector<int8_t> BLOSUM62 = {
// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *
   4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1,-1,-1,-4,  // A 
  -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1,-2, 0,-1,-4,  // R 
  -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 4,-3, 0,-1,-4,  // N 
  -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4,-3, 1,-1,-4,  // D 
   0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-1,-3,-1,-4,  // C 
  -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0,-2, 4,-1,-4,  // Q 
  -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1,-3, 4,-1,-4,  // E 
   0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-4,-2,-1,-4,  // G 
  -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0,-3, 0,-1,-4,  // H 
  -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3, 3,-3,-1,-4,  // I 
  -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4, 3,-3,-1,-4,  // L 
  -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0,-3, 1,-1,-4,  // K 
  -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3, 2,-1,-1,-4,  // M 
  -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3, 0,-3,-1,-4,  // F 
  -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-3,-1,-1,-4,  // P 
   1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0,-2, 0,-1,-4,  // S 
   0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1,-1,-1,-4,  // T 
  -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,-2,-1,-4,  // W 
  -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-1,-2,-1,-4,  // Y 
   0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3, 2,-2,-1,-4,  // V 
  -2,-1, 4, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,-3, 0,-1,-4,  // B 
  -1,-2,-3,-3,-1,-2,-3,-4,-3, 3, 3,-3, 2, 0,-3,-2,-1,-2,-1, 2,-3, 3,-3,-1,-4,  // J 
  -1, 0, 0, 1,-3, 4, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-2,-2,-2, 0,-3, 4,-1,-4,  // Z 
  -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-4,  // X 
  -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,  // * 
};

const swatlib::Matrix<int8_t> NA_MATRIX(16, 16, NA_SUBSTITUTION);
const swatlib::Matrix<int8_t> BLOSUM50_MATRIX(25, 25, BLOSUM50);
const swatlib::Matrix<int8_t> BLOSUM62_MATRIX(25, 25, BLOSUM62);

/**
 * Scale a 1-0 matrix to the given positive and negative values.
 */
template<typename T>
swatlib::Matrix<T> scaleMatrix(swatlib::Matrix<T>, T pos, T neg, T ambig);

Matrix<int8_t> createQueryProfile(const std::vector<uint8_t>&, const swatlib::Matrix<int8_t>&);

enum class Similarity {
    nucleicAcid,  // for na
    blosum50,  // for aa 
    blosum62,  // for aa
    simple,  // simple match mismatch 1/-1
};

Similarity strToSimilarity(std::string);
std::string similarityToStr(Similarity);

template <typename T>
swatlib::Matrix<T> convertMatrix(swatlib::Matrix<int8_t>);

/**
 * Select Matrix while automatically applying match and mismatch values for nucleicAcid similarity table.
 */
swatlib::Matrix<int8_t> selectMatrix(Similarity, int8_t matchValue, int8_t mismatchValue, int8_t ambiguityValue);

}

#endif // SIMILARITY_H