""" Tools for combining images to create videos """

import os
import cv2

def create_video(image_folder_path, video_output_path, frame_rate):
    """
    Combines folder of .png images to create a .mp4 video

    :param image_folder_path: Path to folder containing images.
    :type image_folder_path: str
    :param video_output_path: Destination path for output video
    :type video_output_path: str
    :param frame_rate: Video frames per second
    :type frame_rate: int
    """
    
    images = [img for img in sorted(os.listdir(image_folder_path)) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder_path, images[0]))
    
    height, width, _ = frame.shape
    
    video = cv2.VideoWriter(video_output_path, 0, frame_rate, (width, height))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder_path, image)))
    
    cv2.destroyAllWindows()
    video.release()