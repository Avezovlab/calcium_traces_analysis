# calcium_traces_analysis

Steps to use the code:

1. extract videos from incucyte (uncompressed AVI, best quality)
2. extract the corresponding object masks as b&w tif images (objects being in white).
3. In ImageJ, run the convert_raw_videos.py to create tif stacks from avi files.
4. Run extract_intensity.py to extract active object's locations and calcium traces.
5. Run your analyses (for an example see display_traces.py)
