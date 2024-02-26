from fibers.compose.pipeline_code import tool_box
from fibers.compose.pipeline_code.instruction_runner import InstructionRunner
from fibers.compose.pipeline_text.tree_preprocess import preprocess_text_tree
from fibers.data_loader.markdown_to_tree import markdown_to_tree
from fibers.helper.cache.cache_service import caching
import mentpy as mp

def main():

    import mentpy_tools

    instructions = """
# Lie algebra calculation of a MBQC circuit
- Generate an MBQC circuit of a line cluster state with 5 qubits
- Calculate the Lie algebra of the MBQC circuit
"""

    tree = markdown_to_tree(instructions)

    preprocess_text_tree(tree, fat_limit=300)

    inst_runner = InstructionRunner([mentpy_tools, tool_box], None, external_modules=[("mp", mp)])

    inst_runner.grow_instruction_tree(tree.root)

    caching.save_used()


if __name__ == '__main__':
    main()