@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" goto install

where py >nul 2>nul
if not errorlevel 1 py -3.13 -m venv .venv
if exist ".venv\Scripts\python.exe" goto install

where python >nul 2>nul
if errorlevel 1 goto missing_python
python -m venv .venv
if exist ".venv\Scripts\python.exe" goto install

goto failed

:install
if exist ".venv\hemoflow_dependencies_ready" goto run
echo Installing HemoFlow's Python packages for the first launch...
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto failed
type nul > ".venv\hemoflow_dependencies_ready"

:run
".venv\Scripts\python.exe" main.py %*
if errorlevel 1 goto failed
exit /b 0

:missing_python
echo.
echo HemoFlow could not find Python. Install Python 3.13 from python.org or the Microsoft Store,
echo then open a new Command Prompt and run this launcher again.
pause
exit /b 1

:failed
echo.
echo HemoFlow could not start. The error is shown above.
echo You can also use an existing Python environment:
echo python -m pip install -r requirements.txt
echo python main.py
pause
exit /b 1
