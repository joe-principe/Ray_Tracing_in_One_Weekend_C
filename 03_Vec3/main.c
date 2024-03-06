#include <stdio.h>
#include <stdlib.h>

#include <cglm/cglm.h>

#define IMG_WIDTH 256
#define IMG_HEIGHT 256

int
main(void)
{
    int i, j;
    int ir, ig, ib;
    float r, g, b;

    FILE *img_fp = NULL;
    const char *filename = "image.ppm";

    if ((img_fp = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Error: Could not open %s for writing.\n", filename);
        exit(EXIT_FAILURE);
    }

    fprintf(img_fp, "P3\n%d %d\n255\n", IMG_WIDTH, IMG_HEIGHT);

    for (i = 0; i < IMG_HEIGHT; i++) {
        printf("Scanlines remaining: %d\n", IMG_HEIGHT - i);

        for (j = 0; j < IMG_WIDTH; j++) {
            r = (float)j / (IMG_HEIGHT - 1);
            g = (float)i / (IMG_WIDTH - 1);
            b = 0;

            ir = (int)(255.99 * r);
            ig = (int)(255.99 * g);
            ib = (int)(255.99 * b);

            fprintf(img_fp, "%d %d %d\n", ir, ig, ib);
        }
    }
    
    printf("Done\n");

    fclose(img_fp);

    return 0;
}
/* EOF */
