#!/bin/bash

# Customize your branch and size limit here
REMOTE="origin"
BRANCH="main"
MAXSIZE=$((10 * 1024 * 1024))  # 10 MB

echo "🔍 Checking files to be pushed to $REMOTE/$BRANCH ..."

# Get list of files changed between local branch and remote
FILES=$(git diff --name-only $REMOTE/$BRANCH..HEAD)

if [ -z "$FILES" ]; then
  echo "✅ No files to be pushed."
  exit 0
fi

echo ""
echo "📄 Files to be pushed:"
echo "$FILES"
echo ""

# Check sizes
LARGE_FOUND=0
echo "⚠️ Checking for large files (> $(($MAXSIZE / 1024 / 1024)) MB):"
for FILE in $FILES; do
  # Only check files that exist in the working directory
  if [ -f "$FILE" ]; then
    SIZE=$(stat -c%s "$FILE")
    if [ $SIZE -gt $MAXSIZE ]; then
      echo "🚨 $FILE - $(($SIZE / 1024 / 1024)) MB"
      LARGE_FOUND=1
    fi
  fi
done

if [ $LARGE_FOUND -eq 1 ]; then
  echo ""
  echo "❌ Warning: Some files exceed the size threshold and may cause push issues."
else
  echo "✅ No large files detected."
fi
