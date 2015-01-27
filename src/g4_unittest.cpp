
#include "gtest/gtest.h"
#include "g4.h"



// QGRS Mapper:  45 set to max, gscore = 108
// QGRS Mapper:  30 set to max, gscore = 63
TEST(GScoreTest, 4Tetrad108) {
  G4Candidate c(4, 1, 1, 1);
  cout << "Testing " << c.sequence << endl;
  EXPECT_EQ(c.score_mapper(45), 108);
  EXPECT_EQ(c.score_mapper(30), 63);
  EXPECT_EQ(c.score(), 84);
}

// GGGAGGGAAAAAAAAAAGGGAAGGG
// QGRS Mapper:  45 set to max, gscore = 63
// QGRS Mapper:  30 set to max, gscore = 33
TEST(GScoreTest, 3Tetrad63) {
  G4Candidate c(3, 1, 10, 2);
  cout << "Testing " << c.sequence << endl;
  EXPECT_EQ(c.score_mapper(45), 63);
  EXPECT_EQ(c.score_mapper(30), 33);
  EXPECT_EQ(c.score(), 58);
}

// GGAAGGAAAGGAAAGG
// QGRS Mapper:  45 set to max, gscore = 35
// QGRS Mapper:  30 set to max, gscore = 20
TEST(GScoreTest, 2Tetrad35) {
  G4Candidate c(2, 2, 3, 3);
  cout << "Testing " << c.sequence << endl;
  EXPECT_EQ(c.score_mapper(45), 35);
  EXPECT_EQ(c.score_mapper(30), 20);
  EXPECT_EQ(c.score(), 20);
}

// GGAAAAAAAAAAGGAAAAAGGAAAGG
// QGRS Mapper:  45 set to max, gscore = 29
// QGRS Mapper:  30 set to max, gscore = 14
TEST(GScoreTest, 2Tetrad29) {
  G4Candidate c(2, 10, 5, 3);
  cout << "Testing " << c.sequence << endl;
  EXPECT_EQ(c.score_mapper(45), 29);
  EXPECT_EQ(c.score_mapper(30), 14);
  EXPECT_EQ(c.score(), 16);
}


// GGGGTTGTGGTCAGTGCGGGCCCATGGCCTGGCTGGGCCCGGG
TEST(GScoreTest, 3Tetrad56) {
  G4Candidate c(3, 14, 14, 3);
  cout << "Testing " << c.sequence << endl;
  EXPECT_EQ(c.score(), 56);
}
