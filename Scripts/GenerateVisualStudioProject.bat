@echo off

cd /d %~dp0..

set ARCHITECTURE=x64
set BUILD_DIR=.build

:menu
echo.
echo Select Visual Studio version:
echo   [1] Visual Studio 2022
echo   [2] Visual Studio 2024
echo   [3] Visual Studio 2026
echo   [0] Exit
echo.
set /p CHOICE="Enter choice: "

if "%CHOICE%"=="1" set GENERATOR=Visual Studio 17 2022
if "%CHOICE%"=="2" set GENERATOR=Visual Studio 18 2024
if "%CHOICE%"=="3" set GENERATOR=Visual Studio 18 2026
if "%CHOICE%"=="0" exit /b 0

if not defined GENERATOR (
    echo.
    echo Invalid choice, try again.
    goto menu
)

echo.
echo Generating with: %GENERATOR%
echo.

if not exist %BUILD_DIR% mkdir %BUILD_DIR%
cmake -S . -B %BUILD_DIR% -G "%GENERATOR%" -A %ARCHITECTURE%

pause