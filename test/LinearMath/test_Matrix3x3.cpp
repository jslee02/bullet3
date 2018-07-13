#include <gtest/gtest.h>

GTEST_TEST(LinearMath, Empty)
{
}

int main(int argc, char **argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
