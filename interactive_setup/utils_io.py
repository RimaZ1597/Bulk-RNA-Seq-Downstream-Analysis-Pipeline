"""
utils_io.py
===========
Utility functions for:
- printing headers/sections
- safe_input with defaults
- menu selection
- writing summary reports
Designed following Nextflow training best practices:
clean UI, reproducible config creation, separation of concerns.
"""
def print_header(msg):
    print(f"\n\n⭐ {msg}\n" + "═" * (len(msg) + 4))
def print_section(msg):
    print(f"\n\n➡ {msg}")
    print("─" * (len(msg) + 6))
def safe_input(prompt, default=None):
    if default is not None:
        full_prompt = f"{prompt} [{default}]: "
    else:
        full_prompt = f"{prompt}: "
    value = input(full_prompt).strip()
    
    # Global exit option
    if value.lower() in ['exit', 'quit', 'q']:
        print("\n👋 Exiting setup wizard...")
        exit(0)
    
    return value if value else default
def menu_select(options, prompt="Select option", allow_none=False):
    print(f"\n{prompt}:")
    for i, option in enumerate(options, 1):
        print(f" {i}. {option}")
    if allow_none:
        print(f" {len(options) + 1}. None (skip)")
    while True:
        choice = input("→ Select option: ").strip().lower()
        
        # Global exit option
        if choice in ['exit', 'quit', 'q']:
            print("\n👋 Exiting setup wizard...")
            exit(0)
        
        if choice.isdigit():
            choice_num = int(choice)
            if 1 <= choice_num <= len(options):
                return options[choice_num - 1]
            elif allow_none and choice_num == len(options) + 1:
                return None
        elif choice in ['none', 'skip', 'n'] and allow_none:
            return None
        print("✗ Invalid choice. Try again.")
def confirm(message):
    while True:
        answer = input(f"{message} (y/n): ").strip().lower()
        
        # Global exit option
        if answer in ['exit', 'quit', 'q']:
            print("\n👋 Exiting setup wizard...")
            exit(0)
        
        if answer in ("y", "yes"):
            return True
        if answer in ("n", "no"):
            return False
        print("✗ Enter y/n")
def write_summary_report(path, config):
    with open(path, "w") as f:
        f.write("RNA-seq Analysis Summary\n")
        f.write("=======================\n\n")
        f.write(f"Project: {config['project_name']}\n\n")
        f.write("Experimental Design:\n")
        f.write(f" Primary factor: {config['experimental_design']['primary_factor']}\n")
        f.write(f" Batch factor: {config['experimental_design']['batch_factor']}\n")
        f.write(f" Additional factors: {config['experimental_design']['additional_factors']}\n")
        f.write("Analysis Parameters:\n")
        for k, v in config["analysis_parameters"].items():
            f.write(f" {k}: {v}\n")
        f.write("\nComparisons:\n")
        for c in config["comparisons"]:
            f.write(f" - {c}\n")
    print(f"✓ Summary written: {path}")
