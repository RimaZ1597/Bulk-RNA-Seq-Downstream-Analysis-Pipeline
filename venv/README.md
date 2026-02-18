# Virtual Environment (venv) - Python Dependencies

## Overview
This directory contains Python virtual environments for the BulkRNAseq pipeline. Virtual environments isolate Python dependencies to prevent conflicts between different projects and ensure reproducible execution.

## Virtual Environment Structure

### Main Environment Location
```
~/.venv/                    # Primary virtual environment
├── bin/                    # Executable scripts
│   ├── python              # Python interpreter
│   ├── pip                 # Package installer
│   └── activate            # Activation script
├── lib/                    # Python libraries
│   └── python3.x/
│       └── site-packages/  # Installed packages
├── include/                # Header files
├── pyvenv.cfg             # Environment configuration
└── .gitignore             # Git ignore rules
```

### Project-Specific Environments (if present)
```
project_folder/venv/        # Project-local environment
├── bin/
├── lib/
└── pyvenv.cfg
```

## Required Python Packages

### Core Dependencies (requirements.txt)
```pip
plutobio>=0.1.13           # Pluto.jl integration for interactive analysis
requests>=2.28.0           # HTTP library for API calls
pandas>=1.5.0              # Data manipulation and analysis
pyyaml>=6.0                # YAML configuration file parsing
```

### Additional Dependencies (auto-installed)
- **numpy**: Numerical computing (pandas dependency)
- **urllib3**: HTTP client (requests dependency)
- **certifi**: SSL certificate bundle
- **charset-normalizer**: Unicode text detection
- **idna**: Internationalized domain names

## Environment Management

### Activation/Deactivation
```bash
# Activate virtual environment
source ~/.venv/bin/activate

# Verify activation (should show venv path)
which python
which pip

# Deactivate when done
deactivate
```

### Package Management
```bash
# Install requirements
pip install -r requirements.txt

# Install individual packages
pip install package_name

# List installed packages
pip list

# Show package information
pip show package_name

# Upgrade packages
pip install --upgrade package_name

# Generate requirements file
pip freeze > requirements.txt
```

### Environment Verification
```bash
# Check Python version
python --version

# Verify key packages
python -c "import pandas; print(f'pandas: {pandas.__version__}')"
python -c "import yaml; print(f'PyYAML: {yaml.__version__}')"
python -c "import requests; print(f'requests: {requests.__version__}')"

# Check Pluto integration
python -c "import plutobio; print(f'plutobio: {plutobio.__version__}')"
```

## Pipeline Integration

### Interactive Setup Usage
The virtual environment is primarily used by:
1. **setup.py**: Configuration wizard for analysis parameters
2. **utils_pluto.py**: Pluto.jl notebook integration
3. **explore_pluto.py**: Interactive data exploration
4. **test scripts**: Pipeline validation and testing

### Activation in Scripts
```bash
#!/bin/bash
# Typical pipeline script header
source ~/.venv/bin/activate
python interactive_setup/setup.py
deactivate
```

### Container vs Virtual Environment
- **Development**: Use virtual environment for interactive setup and testing
- **Production**: Use containers (Docker/Singularity) for R analysis scripts
- **Hybrid**: Python components use venv, R components use containers

## Environment Creation and Setup

### Creating New Environment
```bash
# Create new virtual environment
python3 -m venv ~/.venv

# Activate environment
source ~/.venv/bin/activate

# Upgrade pip
pip install --upgrade pip

# Install requirements
pip install -r requirements.txt
```

### Project-Specific Environment
```bash
# Create project-local environment
cd /path/to/project
python3 -m venv venv

# Activate local environment
source venv/bin/activate

# Install requirements
pip install -r requirements.txt
```

### Environment Replication
```bash
# Export current environment
pip freeze > requirements_full.txt

# Create identical environment elsewhere
python3 -m venv new_env
source new_env/bin/activate
pip install -r requirements_full.txt
```

## Troubleshooting

### Common Issues

#### 1. Permission Errors
```bash
# Fix permissions
chmod +x ~/.venv/bin/activate
chmod +x ~/.venv/bin/python
```

#### 2. Package Installation Failures
```bash
# Update pip
pip install --upgrade pip setuptools wheel

# Clear pip cache
pip cache purge

# Install with verbose output
pip install -v package_name
```

#### 3. Import Errors
```bash
# Check if environment is activated
echo $VIRTUAL_ENV

# Verify package installation
pip show package_name

# Check Python path
python -c "import sys; print(sys.path)"
```

#### 4. Environment Corruption
```bash
# Recreate environment
deactivate
rm -rf ~/.venv
python3 -m venv ~/.venv
source ~/.venv/bin/activate
pip install -r requirements.txt
```

### Version Conflicts
```bash
# Check for conflicts
pip check

# Show dependency tree
pip install pipdeptree
pipdeptree

# Force reinstall
pip install --force-reinstall package_name
```

## Best Practices

### Environment Management
1. **Always Activate**: Ensure environment is active before running Python scripts
2. **Regular Updates**: Keep packages updated for security and features
3. **Requirements Tracking**: Maintain accurate requirements.txt file
4. **Environment Testing**: Test environment after major changes

### Development Workflow
```bash
# Start development session
source ~/.venv/bin/activate

# Work on pipeline
cd BulkRNAseq_Pipeline_Limma_Copy
python interactive_setup/setup.py

# Test changes
python -m pytest tests/

# End session
deactivate
```

### Production Deployment
1. **Container Strategy**: Use containers for production R analysis
2. **Version Pinning**: Pin exact package versions in production
3. **Environment Isolation**: Separate environments for dev/test/prod
4. **Backup Requirements**: Keep backup of working requirements.txt

## Security Considerations

### Package Security
```bash
# Check for vulnerabilities
pip install safety
safety check

# Audit packages
pip audit
```

### Environment Isolation
- Never install packages system-wide when using virtual environments
- Use separate environments for different projects
- Regularly update packages to patch security vulnerabilities
- Avoid running untrusted code in shared environments

## Integration with Pipeline Components

### Setup Wizard Integration
```python
# In setup.py
import sys
import os

# Verify virtual environment
if not hasattr(sys, 'real_prefix') and not sys.base_prefix != sys.prefix:
    print("Warning: Virtual environment not detected")
    
# Import pipeline-specific modules
from utils_pluto import load_pluto_data
from utils_io import print_header
```

### Pluto.jl Integration
```python
# PlutoBio package enables Julia integration
import plutobio

# Start Pluto server for interactive analysis
notebook = plutobio.create_notebook("analysis.jl")
notebook.run()
```

## Maintenance Schedule

### Regular Tasks
- **Weekly**: Check for package updates with `pip list --outdated`
- **Monthly**: Update critical packages with `pip install --upgrade`
- **Quarterly**: Recreate environment from clean requirements.txt
- **As Needed**: Add new dependencies and update requirements.txt

### Version Control Integration
```bash
# .gitignore entries for virtual environments
venv/
.venv/
__pycache__/
*.pyc
*.pyo
```

## Monitoring and Logging

### Environment Health Checks
```bash
# Create health check script
cat > check_env.py << 'EOF'
import sys
import pkg_resources

print(f"Python: {sys.version}")
print(f"Virtual env: {sys.prefix}")

required = ['pandas', 'pyyaml', 'requests', 'plutobio']
for pkg in required:
    try:
        dist = pkg_resources.get_distribution(pkg)
        print(f"✓ {pkg}: {dist.version}")
    except pkg_resources.DistributionNotFound:
        print(f"✗ {pkg}: NOT INSTALLED")
EOF

python check_env.py
```

### Usage Tracking
```bash
# Log environment usage
echo "$(date): Environment activated" >> ~/.venv/usage.log
```

## Support and Resources

### Getting Help
1. **Environment Issues**: Recreate from requirements.txt
2. **Package Problems**: Check pip documentation
3. **Integration Issues**: Verify all dependencies are installed
4. **Performance Issues**: Consider using pip-tools for dependency resolution

### Useful Commands Reference
```bash
# Environment info
pip list                    # List all packages
pip show package           # Show package details
pip check                  # Check for broken dependencies
which python              # Show Python executable path
echo $VIRTUAL_ENV         # Show venv path

# Package management
pip install package       # Install package
pip uninstall package     # Remove package
pip freeze               # Export requirements
pip install -r file      # Install from requirements
```

## Version Information
- **Python Version**: 3.8+
- **pip Version**: Latest stable
- **Virtual Environment**: venv (standard library)
- **Last Updated**: December 2025

This virtual environment setup ensures consistent Python dependencies for the interactive components of the BulkRNAseq pipeline while maintaining isolation from system packages.
