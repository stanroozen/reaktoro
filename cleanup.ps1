# Get initial size of build-msvc
Write-Host "Calculating initial size..." -ForegroundColor Yellow
$initialSize = (Get-ChildItem -Path .\build-msvc -Recurse -Force -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum
Write-Host "Initial build-msvc size: $([math]::Round($initialSize / 1GB, 3)) GB"

# Define file patterns to delete
$patterns = @('*.ilk', '*.ipdb', '*.iobj', '*.pdb')

# Find and delete files
Write-Host "`nDeleting build artifacts..." -ForegroundColor Yellow
$totalDeleted = 0
$fileCount = 0

foreach ($pattern in $patterns) {
    $files = Get-ChildItem -Path .\build-msvc -Filter $pattern -Recurse -Force -ErrorAction SilentlyContinue
    if ($files) {
        Write-Host "`nFiles matching $pattern :"
        foreach ($file in $files) {
            $sizeMB = [math]::Round($file.Length / 1MB, 2)
            Write-Host "  $($file.FullName.Replace($PWD, '.')) - $sizeMB MB"
            Remove-Item -Path $file.FullName -Force -ErrorAction SilentlyContinue
            $totalDeleted += $file.Length
            $fileCount++
        }
    }
}

# Get final size of build-msvc
Write-Host "`nCalculating final size..." -ForegroundColor Yellow
$finalSize = (Get-ChildItem -Path .\build-msvc -Recurse -Force -ErrorAction SilentlyContinue | Measure-Object -Property Length -Sum).Sum

Write-Host "`n========== CLEANUP COMPLETE ==========" -ForegroundColor Green
Write-Host "Files deleted: $fileCount"
Write-Host "Total space reclaimed: $([math]::Round($totalDeleted / 1GB, 3)) GB ($([math]::Round($totalDeleted / 1MB, 2)) MB)"
Write-Host "build-msvc size before: $([math]::Round($initialSize / 1GB, 3)) GB"
Write-Host "build-msvc size after: $([math]::Round($finalSize / 1GB, 3)) GB"
Write-Host "build-msvc reduced by: $([math]::Round(($initialSize - $finalSize) / 1GB, 3)) GB"

# Get remaining free space on C:
$driveInfo = Get-PSDrive -Name C
$freeSpace = $driveInfo.Free
Write-Host "`nRemaining free space on C: : $([math]::Round($freeSpace / 1GB, 2)) GB"
Write-Host "========================================" -ForegroundColor Green
