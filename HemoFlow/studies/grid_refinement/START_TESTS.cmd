@echo off
setlocal
pushd "%~dp0"
set "HEMO_PY=.venv\Scripts\python.exe"
if not exist "%HEMO_PY%" set "HEMO_PY=..\.venv\Scripts\python.exe"
if not exist "%HEMO_PY%" set "HEMO_PY=..\..\.venv\Scripts\python.exe"
if not exist "%HEMO_PY%" set "HEMO_PY=python"
"%HEMO_PY%" -c "import numpy, numba" >nul 2>&1
if errorlevel 1 (
  echo Python with NumPy and Numba was not found.
  echo Use your HemoFlow environment, or install requirements.txt into an environment.
  pause
  popd
  exit /b 1
)
"%HEMO_PY%" run_grid.py %*
set "HEMO_EXIT=%ERRORLEVEL%"
echo.
echo HemoFlow runner exit code: %HEMO_EXIT%
pause
popd
exit /b %HEMO_EXIT%
