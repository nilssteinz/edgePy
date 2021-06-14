import os
import unittest

print(os.path.split(os.getcwd()))
class test_edgePy(unittest.TestCase):
    def setUp(self) -> None:
        print(os.path.basename(os.getcwd()))
        pass

    def tearDown(self) -> None:
        pass

    def test_stuff(self):
        self.assertEqual(1,2)

if __name__ == "__main__":
    unittest.main()
