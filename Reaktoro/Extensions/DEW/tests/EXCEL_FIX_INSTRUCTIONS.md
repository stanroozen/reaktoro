# How to Fix the Excel VBA Export Code for depsdrho Truth Tables

## The Problem

The `calculate_depsdrho` function **modifies the `rho` parameter in-place**, converting it from g/cm³ to mol/cm³. When the export code passes `rho` to this function, the variable gets changed from (e.g.) `0.1` to `0.0055508`.

Sequence of events:
1. `rho = 0.1` (from RHO_LIST)
2. `v = calculate_depsdrho(rho, T, eqe)` - **Inside this function, `rho` is multiplied by 0.055508**
3. After function returns: `rho = 0.0055508` (modified!)
4. `DotNum(rho, ...)` writes the **modified** value to CSV

This creates a mismatch where:
- CSV column 3 shows: `0.0055508` (modified value after function call)
- But the calculation used the original: `0.1` g/cm³

## The Fix

Save the original `rho` value in a separate variable **before** calling `calculate_depsdrho()`, then use that saved value when writing to CSV.

### Current Buggy Code:

```vba
Private Sub Export_depsdrho_All()
    EnsureGrids
    On Error Resume Next

    Dim rows() As String, n As Long
    Dim Ti As Long, eqe As Integer, r As Long
    Dim T As Double, rho As Double, v As Double
    Dim path As String, tag As String
    Dim RHO_LIST As Variant
    Dim vStr As String

    RHO_LIST = Array(0.1, 0.5, 0.9)

    For eqe = 1 To 4
        tag = Choose(eqe, "JN1991", "Franck1990", "Fernandez1997", "Power")
        path = OutPath("truth_depsdrho_" & tag & ".csv")
        InitBuf rows, 128: n = 0
        AppendRow rows, n, "T_C,epsEq_id,rho_g_per_cm3,depsdrho_cm3_per_g"

        For Ti = LBound(T_LIST) To UBound(T_LIST)
            T = T_LIST(Ti)
            If ValidDepsdrho(T, eqe) Then
                For r = LBound(RHO_LIST) To UBound(RHO_LIST)
                    rho = CDbl(RHO_LIST(r))  ' rho = 0.1
                    v = 0: vStr = ""

                    Err.Clear
                    v = calculate_depsdrho(rho, T, eqe)  ' <-- rho becomes 0.0055508!
                    If Err.Number = 0 Then
                        vStr = DotNum(v, "0.000000000000E+00")
                    Else
                        Err.Clear
                    End If

                    ' BUG: rho has been modified by calculate_depsdrho!
                    AppendRow rows, n, _
                        DotNum(T, "0.00") & "," & CStr(eqe) & "," & _
                        DotNum(rho, "0.000000000000") & "," & vStr
                        ' ^^^^ writes 0.0055508 instead of 0.1
                Next r
            End If
        Next Ti

        FlushCSV path, rows, n
    Next eqe

    On Error GoTo 0
End Sub
```

### Corrected Code:

```vba
Private Sub Export_depsdrho_All()
    EnsureGrids
    On Error Resume Next

    Dim rows() As String, n As Long
    Dim Ti As Long, eqe As Integer, r As Long
    Dim T As Double, rho As Double, v As Double
    Dim path As String, tag As String
    Dim RHO_LIST As Variant
    Dim vStr As String
    Dim rho_original As Double  ' <-- ADD THIS: Save original rho value

    RHO_LIST = Array(0.1, 0.5, 0.9)

    For eqe = 1 To 4
        tag = Choose(eqe, "JN1991", "Franck1990", "Fernandez1997", "Power")
        path = OutPath("truth_depsdrho_" & tag & ".csv")
        InitBuf rows, 128: n = 0
        AppendRow rows, n, "T_C,epsEq_id,rho_g_per_cm3,depsdrho_cm3_per_g"

        For Ti = LBound(T_LIST) To UBound(T_LIST)
            T = T_LIST(Ti)
            If ValidDepsdrho(T, eqe) Then
                For r = LBound(RHO_LIST) To UBound(RHO_LIST)
                    rho = CDbl(RHO_LIST(r))
                    rho_original = rho  ' <-- SAVE original value before function call
                    v = 0: vStr = ""

                    Err.Clear
                    v = calculate_depsdrho(rho, T, eqe)  ' rho gets modified here
                    If Err.Number = 0 Then
                        vStr = DotNum(v, "0.000000000000E+00")
                    Else
                        Err.Clear
                    End If

                    ' FIXED: Write the original rho value, not the modified one
                    AppendRow rows, n, _
                        DotNum(T, "0.00") & "," & CStr(eqe) & "," & _
                        DotNum(rho_original, "0.000000000000") & "," & vStr
                        ' ^^^^^^^^^^^^ Use saved original value
                Next r
            End If
        Next Ti

        FlushCSV path, rows, n
    Next eqe

    On Error GoTo 0
End Sub
```

## Expected Output After Fix

### Before (buggy):
```
T_C,epsEq_id,rho_g_per_cm3,depsdrho_cm3_per_g
25.00,2,0.005550800000,3.988988482479E+01
25.00,2,0.027754000000,1.119881396419E+02
25.00,2,0.049957200000,1.365857052352E+03
```

### After (correct):
```
T_C,epsEq_id,rho_g_per_cm3,depsdrho_cm3_per_g
25.00,2,0.100000000000,3.988988482479E+01
25.00,2,0.500000000000,1.119881396419E+02
25.00,2,0.900000000000,1.365857052352E+03
```

## Why This Matters

The CSV header says `rho_g_per_cm3` which means the values should be in **g/cm³**, but the buggy code writes values that are in strange mixed units (g/cm³ × mol/g = mol/cm³).

After the fix:
- Column 3 will contain actual density in g/cm³ (0.1, 0.5, 0.9)
- Column 4 will contain the depsdrho calculated at those densities
- The C++ code will correctly read and compare against these physically meaningful values

## Steps to Regenerate Truth Tables

1. Find and fix the `Export_depsdrho_All()` function in your Excel VBA code
2. Remove any `* 0.055508` multiplication before writing density to CSV
3. Run `Export_All_TruthCSVs()` to regenerate all truth tables
4. Copy the new `truth_depsdrho_*.csv` files to the C++ test directory
5. Rebuild and run the C++ tests - they should now pass!
