
#include "gtest/gtest.h"
#include "g4.h"

TEST(GScoreTest, 4Tetrad84) {
  G4Candidate c(4, 1, 1, 1);
  EXPECT_EQ(c.score(), 84);
}

TEST(GScoreTest, 3Tetrad58) {
  G4Candidate c(3, 1, 10, 2);
  EXPECT_EQ(c.score(), 58);
}

TEST(GScoreTest, 2Tetrad20) {
  G4Candidate c(2, 2, 3, 3);
  EXPECT_EQ(c.score(), 20);
}
