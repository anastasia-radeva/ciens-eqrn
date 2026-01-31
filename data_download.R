options(timeout = 24*3600)  # 24h

urls <- c(
  "https://radar.kit.edu/radar-backend/archives/zzfEJPxDILXwSNPH/versions/1/content",
  "https://radar.kit.edu/radar-backend/archives/JKALdQqqLIjGUOBC/versions/1/content"
)

dest_00 <- "~/CIENS_Data/10.35097-zzfEJPxDILXwSNPH.tar"
dest_12 <- "~/CIENS_Data/10.35097-JKALdQqqLIjGUOBC.tar"

dir.create(dirname(dest_12), recursive = TRUE, showWarnings = FALSE)
download.file(
  urls[2], destfile = dest_12, mode = "wb",
  method = "curl",
  extra = "-L -C - --retry 20 --retry-delay 5 --retry-connrefused --speed-limit 100000 --speed-time 60"
)

dir.create(dirname(dest_00), recursive = TRUE, showWarnings = FALSE)
download.file(
  urls[1], destfile = dest_00, mode = "wb",
  method = "curl",
  extra = "-L -C - --retry 20 --retry-delay 5 --retry-connrefused --speed-limit 100000 --speed-time 60"
)

# -------------------------------------------------------------------
# Extract the tar files
# -------------------------------------------------------------------
extract_00_dir <- "~/CIENS_Data/10.35097-zzfEJPxDILXwSNPH"
extract_12_dir <- "~/CIENS_Data/10.35097-JKALdQqqLIjGUOBC"

dir.create(extract_00_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_12_dir, recursive = TRUE, showWarnings = FALSE)

utils::untar(dest_00, exdir = extract_00_dir)
utils::untar(dest_12, exdir = extract_12_dir)

# -------------------------------------------------------------------
# Copy run0, run12, observations to /CIENS_Data/data/dataset/
# -------------------------------------------------------------------
target_root <- "~/CIENS_Data/data/dataset"
dir.create(target_root, recursive = TRUE, showWarnings = FALSE)

copy_if_exists <- function(src, dst) {
  if (dir.exists(src)) {
    # If destination exists, remove it first to avoid nested copies
    if (dir.exists(dst)) unlink(dst, recursive = TRUE, force = TRUE)
    ok <- file.copy(src, dst, recursive = TRUE)
    if (!ok) warning("Failed to copy: ", src, " -> ", dst)
  } else {
    warning("Source folder not found (skipping): ", src)
  }
}

# Sources
srcs <- list(
  list(root = extract_00_dir, name = "run0"),
  list(root = extract_00_dir, name = "observations"),
  list(root = extract_12_dir, name = "run12")
  )

move_if_exists <- function(src, dst) {
  if (!dir.exists(src)) {
    warning("Source folder not found (skipping): ", src)
    return(invisible(FALSE))
  }
  
  # If destination exists, remove it first to avoid nesting/partial merges
  if (dir.exists(dst)) unlink(dst, recursive = TRUE, force = TRUE)
  
  # Try fast move first (works if same filesystem)
  ok <- file.rename(src, dst)
  
  # Fallback: copy + delete (works across filesystems)
  if (!ok) {
    ok_copy <- file.copy(src, dst, recursive = TRUE)
    if (!ok_copy) {
      warning("Failed to move (copy step failed): ", src, " -> ", dst)
      return(invisible(FALSE))
    }
    unlink(src, recursive = TRUE, force = TRUE)
    ok <- TRUE
  }
  
  invisible(ok)
}

# Move
for (s in srcs) {
  src_path <- file.path(s$root, s$name)
  dst_path <- file.path(target_root, s$name)
  move_if_exists(src_path, dst_path)
}

message("Done. Moved folders into:\n  ", target_root)
