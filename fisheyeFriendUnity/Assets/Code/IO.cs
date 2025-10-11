using System.IO;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.EventSystems;
using SFB;

namespace FisheyeFriend {
    public class IO : MonoBehaviour {


        // Assign in the Inspector: Button for opening file dialog
        public Button openImageButton; // "Open Image" button

        // Assign in the Inspector: Button for saving the original image
        public Button saveImageButton; // "Save" button

        // Assign in the Inspector: RawImage for displaying the resized image (left)
        public RawImage leftImage; // Left RawImage

        // Assign in the Inspector: RawImage for displaying the resized image (right)
        public RawImage rightImage; // Right RawImage

        // Assign in the Inspector: ScrollRect containing both images
        public ScrollRect scrollRect; // ScrollRect for side-by-side images

        // Store the original loaded image
        private Texture2D originalTexture;

        // Store the resized image
        private Texture2D resizedTexture;

        void Start () {
            openImageButton.onClick.AddListener(OnOpenImageClicked);
            saveImageButton.onClick.AddListener(OnSaveImageClicked);
            saveImageButton.interactable = false; // Disable until an image is loaded
        }

        void OnOpenImageClicked () {
            Debug.Log("Open Image button clicked.");
            // Use a file picker plugin or a custom file dialog for Unity.
            // For demonstration, this is a placeholder path.
            string path = ShowFilePicker(); // Implement this with a plugin or native dialog

            if (!string.IsNullOrEmpty(path) && File.Exists(path)) {
                byte[] fileData = File.ReadAllBytes(path);
                originalTexture = new Texture2D(2, 2);
                if (originalTexture.LoadImage(fileData)) {
                    resizedTexture = ResizeTextureMaintainingAspect(originalTexture, 900, 400);

                    // Display the resized image on both sides
                    leftImage.texture = resizedTexture;
                    rightImage.texture = Instantiate(resizedTexture);

                    // Adjust the size of the RawImages to fit the resized texture
                    Vector2 size = new Vector2(resizedTexture.width, resizedTexture.height);
                    leftImage.rectTransform.sizeDelta = size;
                    rightImage.rectTransform.sizeDelta = size;

                    // Enable the save button
                    saveImageButton.interactable = true;
                }
            }
        }

        void OnSaveImageClicked () {
            // Placeholder for save logic
            SaveOriginalImage();
        }

        // Resize the texture maintaining aspect ratio
        Texture2D ResizeTextureMaintainingAspect (Texture2D source, int maxWidth, int maxHeight) {
            float ratio = Mathf.Min((float)maxWidth / source.width, (float)maxHeight / source.height, 1f);
            int newWidth = Mathf.RoundToInt(source.width * ratio);
            int newHeight = Mathf.RoundToInt(source.height * ratio);

            Texture2D result = new Texture2D(newWidth, newHeight, source.format, false);

            for (int y = 0; y < newHeight; y++) {
                for (int x = 0; x < newWidth; x++) {
                    float u = x / (float)(newWidth - 1);
                    float v = y / (float)(newHeight - 1);
                    Color color = source.GetPixelBilinear(u, v);
                    result.SetPixel(x, y, color);
                }
            }
            result.Apply();
            return result;
        }

        // Placeholder for file picker (implement with a plugin or native dialog)
        string ShowFilePicker () {
            var paths = StandaloneFileBrowser.OpenFilePanel("Open File", "", "", false);
            return paths[0];
        }

        // Placeholder for save logic
        void SaveOriginalImage () {
            // TODO: Implement save dialog and write originalTexture to file
            Debug.Log("SaveOriginalImage called.");
        }


    }
}