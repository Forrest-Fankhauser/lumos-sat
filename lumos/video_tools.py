import os
import cv2

def create_video(image_folder, video_name, frame_rate = 15):
    
    images = [img for img in sorted(os.listdir(image_folder)) if img.endswith(".png")]
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    
    height, width, layers = frame.shape
    
    video = cv2.VideoWriter(video_name, 0, frame_rate, (width, height))
    
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder, image)))
    
    cv2.destroyAllWindows()
    video.release()