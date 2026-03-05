@echo off
echo ========================================
echo QuasiLocal Toolkit - Run Complete Example
echo ========================================
echo.

cd /d C:\Users\PC\Downloads\quasilocal-toolkit\build

echo [1/4] Setting up x64 environment...
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" > nul 2>&1
echo [OK] Environment ready
echo.

echo [2/4] Checking example file...
if not exist simple_example.cpp (
    echo Creating simple_example.cpp...
    
    echo #include ^<iostream^> > simple_example.cpp
    echo #include ^<vector^> >> simple_example.cpp
    echo #include ^<cmath^> >> simple_example.cpp
    echo. >> simple_example.cpp
    echo // Forward declarations >> simple_example.cpp
    echo namespace qlt { >> simple_example.cpp
    echo     std::vector^<double^> compute_axisymmetric_j( >> simple_example.cpp
    echo         const std::vector^<double^>^& alpha, >> simple_example.cpp
    echo         const std::vector^<double^>^& f, >> simple_example.cpp
    echo         const std::vector^<double^>^& r); >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     double compute_axisymmetric_surgical_flux( >> simple_example.cpp
    echo         double alpha1, double r1, double f1, >> simple_example.cpp
    echo         double alpha2, double r2, double f2); >> simple_example.cpp
    echo } >> simple_example.cpp
    echo. >> simple_example.cpp
    echo int main() { >> simple_example.cpp
    echo     std::cout ^<^< "=== QuasiLocal Toolkit Example ===\n\n"; >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     // Test data >> simple_example.cpp
    echo     std::vector^<double^> r = {1.0, 2.0, 3.0, 4.0, 5.0}; >> simple_example.cpp
    echo     std::vector^<double^> alpha = {0.9, 0.95, 0.98, 0.99, 1.0}; >> simple_example.cpp
    echo     std::vector^<double^> f = {0.1, 0.2, 0.3, 0.4, 0.5}; >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     std::cout ^<^< "Input data:\n"; >> simple_example.cpp
    echo     for (size_t i = 0; i ^< r.size(); i++) { >> simple_example.cpp
    echo         std::cout ^<^< "  r = " ^<^< r[i] >> simple_example.cpp
    echo                   ^<^< ", alpha = " ^<^< alpha[i] >> simple_example.cpp
    echo                   ^<^< ", f = " ^<^< f[i] ^<^< "\n"; >> simple_example.cpp
    echo     } >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     // Compute J(r) >> simple_example.cpp
    echo     std::vector^<double^> J = qlt::compute_axisymmetric_j(alpha, f, r); >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     std::cout ^<^< "\nResults J(r):\n"; >> simple_example.cpp
    echo     for (size_t i = 0; i ^< J.size(); i++) { >> simple_example.cpp
    echo         std::cout ^<^< "  J(" ^<^< r[i] ^<^< ") = " ^<^< J[i] ^<^< "\n"; >> simple_example.cpp
    echo     } >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     // Compute surgical flux between r=2 and r=4 >> simple_example.cpp
    echo     double flux = qlt::compute_axisymmetric_surgical_flux( >> simple_example.cpp
    echo         alpha[1], r[1], f[1], >> simple_example.cpp
    echo         alpha[3], r[3], f[3]); >> simple_example.cpp
    echo     std::cout ^<^< "\nSurgical flux between r=2 and r=4: " ^<^< flux ^<^< "\n"; >> simple_example.cpp
    echo. >> simple_example.cpp
    echo     std::cout ^<^< "\n=== Example Complete ===\n"; >> simple_example.cpp
    echo     return 0; >> simple_example.cpp
    echo } >> simple_example.cpp
    
    echo [OK] Example file created
) else (
    echo [OK] Example file exists
)
echo.

echo [3/4] Compiling...
cl /EHsc /MD simple_example.cpp src\core\Release\qlt_core.lib > compile_log.txt 2>&1
if %errorlevel% equ 0 (
    echo [OK] Compilation successful
    type compile_log.txt | findstr /v "Microsoft" | findstr /v "Copyright" | findstr /v "LINK"
) else (
    echo [FAILED] Compilation failed
    echo.
    echo === COMPILATION ERRORS ===
    type compile_log.txt
    echo ===========================
    goto :error
)
echo.

echo [4/4] Running example...
echo.
echo ----------------------------------------
simple_example.exe
echo ----------------------------------------
echo.

echo ========================================
echo [SUCCESS] Example completed successfully!
echo ========================================
goto :end

:error
echo ========================================
echo [FAILED] Example failed!
echo ========================================

:end
pause