# Team Git Workflow 

This is the simple, high-level workflow we follow. It assumes you’re using GitHub Desktop.

## Core Rules
- `main` is always stable and runnable.
- Do not push directly to `main`.
- Every change goes through a Pull Request (PR).

## Branches
Work happens on short-lived branches:
- `feat/<short-description>`
- `fix/<short-description>`
- `refactor/<short-description>`

Examples: `feat/add-csv-import`, `fix/missing-metric`, `refactor/rename-eval`

## Everyday Workflow (GitHub Desktop)
1. **Sync first**  
   In GitHub Desktop, click `Fetch origin`, then `Pull` if it’s available.
2. **Create a branch**  
   `Current branch` → `New branch…` → name it using the patterns above.
3. **Make your changes**  
   Work normally in your editor.
4. **Commit locally**  
   In GitHub Desktop, write a clear summary and click `Commit`.
5. **Publish your branch**  
   Click `Publish branch`.
6. **Open a PR**  
   Click `Create Pull Request` and fill in a short description of what changed and why.
7. **Respond to feedback**  
   If you get comments, update your branch and commit again. The PR updates automatically.
8. **Merge**  
   After approval and checks pass, use **Squash and merge**.
9. **Delete your branch**  
   After merging, delete the branch in GitHub (and in GitHub Desktop when prompted).

## Pull Request Expectations
- Small and focused (one concern per PR).
- Refactors should be mechanical only (renames/moves, no behavior changes).
- If you need refactor + behavior changes, split into two PRs.

## Staying Up-to-Date
Before you open or merge a PR:
- Use GitHub Desktop’s `Fetch origin` and `Pull` to bring in the latest `main`.
- If your branch is behind, GitHub Desktop will offer `Update from main` — use that before asking for review.
