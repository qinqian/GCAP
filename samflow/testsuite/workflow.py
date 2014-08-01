import unittest
from samflow.workflow import Workflow, attach_back, attach_front
from samflow.command import ShellCommand
from samflow.command_on_jinja import JinjaFormatShellCommand as JinShCommand

class WorkflowTestSuite(unittest.TestCase):
    def create_tree(self):
        main_workflow = Workflow("main")
        sub_workflow = Workflow("sub")
        attach_back(sub_workflow, ShellCommand('echo "subtree started"'))
        attach_back(sub_workflow, JinShCommand('touch {{ output|join(" ") }}', output=["f1","f2"]))
        attach_back(sub_workflow, JinShCommand('rm {{ input|join(" ") }}', input=["f1", "f2"]))
        attach_back(sub_workflow, ShellCommand('echo "subtree ended"'))
        attach_back(main_workflow, sub_workflow)
        return main_workflow

    def test_workflow_invoke_success(self):
        tree = self.create_tree()
        self.assertTrue(tree.invoke())
    def test_workflow_attach_later_invoke_success(self):

        tree = self.create_tree()
        attach_front(tree, ShellCommand("touch {output}", output="outer_f1"))
        attach_back(tree, ShellCommand("rm {input}", input="outer_f1"))
        attach_front(tree, ShellCommand('echo "{0} decorator started {0}"'.format("="*10)))
        attach_back(tree, JinShCommand('echo "{0} decorator ended {0}"'.format("="*10)))
        self.assertTrue(tree.invoke())

    def test_workflow_failed_halfway_missing_output(self):
        tree = self.create_tree()
        attach_back(tree, ShellCommand("echo NothingHappened", output="strange_file"))
        self.sleep_after(tree)
        self.assertFalse(tree.invoke())

        tree = self.create_tree()
        attach_back(tree, ShellCommand("echo NothingHappened", output=["s1", "s2"]))
        self.sleep_after(tree)
        self.assertFalse(tree.invoke())

        tree = self.create_tree()
        attach_back(tree, ShellCommand("echo NothingHappened", output={"s1":"s1", "s2":"s2"}))
        self.sleep_after(tree)
        self.assertFalse(tree.invoke())


    def test_workflow_failed_halfway_missing_input(self):
        tree = self.create_tree()
        attach_front(tree, ShellCommand("touch s1", output="s1"))
        attach_back(tree, ShellCommand("rm s1"))
        attach_back(tree, ShellCommand("cat s1", input="s1"))
        self.assertFalse(tree.invoke())

    def sleep_before(self, tree):
        attach_front(tree, ShellCommand("sleep 100"))

    def sleep_after(self, tree):
        attach_back(tree, ShellCommand("sleep 100"))

    def test_workflow_failed_from_start_dangling(self):
        tree = self.create_tree()
        self.sleep_before(tree)
        attach_back(tree, ShellCommand("echo NeedWrongFile", input="wrong_file"))
        self.assertFalse(tree.invoke())

        tree = self.create_tree()
        self.sleep_before(tree)
        attach_front(tree, ShellCommand("echo NeedWrongFileList", input=["wrong_file1", "wrong_file2", "wrong_file3"]))
        self.assertFalse(tree.invoke())

        tree = self.create_tree()
        self.sleep_before(tree)
        attach_back(tree, ShellCommand("echo ala", input={"f1": "f1", "f2": "f2", "f3": "f3"}))
        self.assertFalse(tree.invoke())

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(WorkflowTestSuite))
    return suite
