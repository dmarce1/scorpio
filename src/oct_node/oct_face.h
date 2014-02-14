#ifndef OCT_FACE_H_
#define OCT_FACE_H_

#define OCT_NSIB 6

typedef int OctFace;

#define XL 0
#define XU 1
#define YL 2
#define YU 3
#define ZL 4
#define ZU 5

OctFace invert( OctFace );

#endif /* OCT_FACE_H_ */
