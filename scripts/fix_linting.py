#!/usr/bin/env python3
"""
Script to gradually fix linting issues in the codebase.
This script addresses common, automatable linting violations.
"""

import re
import subprocess
import sys
from pathlib import Path


def fix_trailing_whitespace_and_blank_lines(file_path):
    """Fix trailing whitespace and whitespace-only blank lines."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Remove trailing whitespace from all lines
        content = re.sub(r'[ \t]+$', '', content, flags=re.MULTILINE)
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"Error fixing {file_path}: {e}")
        return False


def fix_f_string_placeholders(file_path):
    """Fix f-strings that are missing placeholders."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original_content = content
        
        # Simple fix: convert f"string" to "string" if no placeholders
        # This is a basic fix - manual review may be needed for complex cases
        lines = content.split('\n')
        fixed_lines = []
        
        for line in lines:
            # Look for f-strings without {} placeholders
            if 'f"' in line and '{' not in line:
                line = line.replace('f"', '"')
            if "f'" in line and '{' not in line:
                line = line.replace("f'", "'")
            fixed_lines.append(line)
        
        content = '\n'.join(fixed_lines)
        
        if content != original_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"Fixed f-string placeholders in {file_path}")
        
        return True
    except Exception as e:
        print(f"Error fixing f-strings in {file_path}: {e}")
        return False


def remove_unused_imports(file_path):
    """Remove obviously unused imports (basic cases only)."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        lines = content.split('\n')
        fixed_lines = []
        
        for line in lines:
            # Remove specific unused imports that we know are safe to remove
            if (line.strip().startswith('from typing import Dict') and 
                'Dict' not in content.replace(line, '')):
                print(f"Removing unused import: {line.strip()} from {file_path}")
                continue
            elif (line.strip().startswith('from typing import Tuple') and 
                  'Tuple' not in content.replace(line, '')):
                print(f"Removing unused import: {line.strip()} from {file_path}")
                continue
            elif line.strip() == 'import os' and 'os.' not in content.replace(line, ''):
                print(f"Removing unused import: {line.strip()} from {file_path}")
                continue
            
            fixed_lines.append(line)
        
        content = '\n'.join(fixed_lines)
        
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"Error removing unused imports in {file_path}: {e}")
        return False


def fix_python_file(file_path):
    """Apply all automated fixes to a Python file."""
    print(f"Fixing: {file_path}")
    
    # Fix whitespace issues
    fix_trailing_whitespace_and_blank_lines(file_path)
    
    # Fix f-string issues
    fix_f_string_placeholders(file_path)
    
    # Remove some unused imports
    remove_unused_imports(file_path)


def main():
    """Main function to fix linting issues."""
    project_root = Path(__file__).parent.parent
    
    # Directories to process
    directories = [
        project_root / 'src',
        project_root / 'tests'
    ]
    
    print("🔧 Starting automated linting fixes...")
    
    # Process all Python files
    for directory in directories:
        if directory.exists():
            for py_file in directory.rglob('*.py'):
                fix_python_file(py_file)
    
    # Fix specific files in root
    for file_name in ['run_tests.py', 'setup.py']:
        file_path = project_root / file_name
        if file_path.exists():
            fix_python_file(file_path)
    
    # Run black to ensure consistent formatting
    print("\n🎨 Running black for consistent formatting...")
    try:
        subprocess.run([
            sys.executable, '-m', 'black',
            str(project_root / 'src'),
            str(project_root / 'tests'),
            str(project_root / 'run_tests.py'),
            str(project_root / 'setup.py')
        ], check=True, cwd=project_root)
        print("✅ Black formatting completed")
    except subprocess.CalledProcessError as e:
        print(f"❌ Black formatting failed: {e}")
    
    print("\n📊 Generating updated linting report...")
    try:
        result = subprocess.run([
            sys.executable, '-m', 'flake8',
            str(project_root / 'src'),
            str(project_root / 'tests'),
            '--count', '--statistics'
        ], capture_output=True, text=True, cwd=project_root)
        
        print("Current linting status:")
        print(result.stdout)
        if result.stderr:
            print("Errors:")
            print(result.stderr)
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Linting check failed: {e}")
    
    print("\n✅ Automated fixes completed!")
    print("Manual review recommended for complex issues.")


if __name__ == '__main__':
    main()
