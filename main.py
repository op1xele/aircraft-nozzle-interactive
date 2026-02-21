"""
Aircraft Nozzle Interactive - Main Entry Point
Advanced Nozzle Design Suite
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gui.main_window import MainWindow


def main():
    print("=" * 60)
    print("Aircraft Nozzle Interactive v2.0")
    print("Pickle Brothers Engineering Team")
    print("=" * 60)
    print()
    
    try:
        app = MainWindow()
        app.run()
    
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())