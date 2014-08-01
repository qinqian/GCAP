import unittest
from samflow.command import ShellCommand

class ShellCommandTestCase(unittest.TestCase):

    def touch_files(self, *files):
        self.assertTrue(ShellCommand("touch " + " ".join(files), output=files).invoke())

    def delete_files(self, *files):
        self.assertTrue(ShellCommand("rm " + " ".join(files), input=files).invoke())

    def test_invoke_success_status(self):
        self.assertTrue(ShellCommand("echo 'run successfully' > /dev/null").invoke())
        self.assertTrue(ShellCommand("exit 0").invoke())

    def test_invoke_fail_status(self):
        self.assertFalse((ShellCommand(template="exit 1")).invoke())

    def test_invoke_collect_output(self):
        echo_cmd = ShellCommand("echo test_collect").set_stdout_collecting()
        self.assertTrue(echo_cmd.invoke())
        self.assertEqual(echo_cmd.result, "test_collect\n")

    def test_invoke_non_exist_input(self):
        non_exist_input_cmd = ShellCommand("cat < {input}", input="non_exist_file")
        self.assertFalse(non_exist_input_cmd.invoke())

    def test_invoke_non_exist_output(self):
        non_exist_output_cmd = ShellCommand("echo tempfile3", output="tempfile3")
        self.assertFalse(non_exist_output_cmd.invoke())


    def test_invoke_exist_input(self):
        self.touch_files("temp_file", "temp_file.1")
        self.delete_files("temp_file", "temp_file.1")

    def test_invoke_dangling_tool(self):
        dangling_tool_command = ShellCommand("{tool} fun", tool="wolfyp")
        self.assertFalse(dangling_tool_command.invoke())

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(ShellCommandTestCase))
    return suite


