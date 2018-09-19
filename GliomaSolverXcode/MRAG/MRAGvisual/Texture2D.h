#pragma once


struct Texture2D
{
	GLuint handle;
	GLint level;
    GLint internalformat;
    GLsizei size[2];
    GLint border;
  //  GLenum format;
    GLenum type;
  //  const GLvoid *pixels;
	GLint filter;

	unsigned int MinP2G(unsigned int val)
	{
		double log_2 = log ((double)val)/log(2.0);
		int shift = (unsigned int)(log_2+0.9999);
		//printf("val=%d res %d \n",val, 1<<shift);
		
		return 1<<shift;
	}
	
	Texture2D( GLsizei width, GLsizei height, GLint _internalformat, GLint _filter, bool bKeepNonPowerOfTwo = false): 
		handle(0),
		level(0),
		internalformat(_internalformat),
		border(0),
		type(0),
		filter(_filter)
	{
		size[0] = bKeepNonPowerOfTwo ? width : MinP2G(width);
		size[1] = bKeepNonPowerOfTwo ? height: MinP2G(height);

		glGenTextures(1,&handle);
		glBindTexture(GL_TEXTURE_2D, handle);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, _filter);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, _filter);
		
		unsigned char * data = NULL;
		
		glTexImage2D(GL_TEXTURE_2D, 0, internalformat, size[0], size[1], 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
		glBindTexture(GL_TEXTURE_2D,0);
	}

	void ShowContent(bool bWriteDepth=false) const 
	{
		glBindTexture(GL_TEXTURE_2D, handle);

		glPushMatrix();
		glLoadIdentity();
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		glOrtho(-1,1,-1,1,0,1);
		glMatrixMode(GL_MODELVIEW);

		glBegin(GL_QUADS);
		{
			glTexCoord2f(0,0);
			glVertex2f(-1,-1); 
				glTexCoord2f(1,0); 
			glVertex2f(1,-1); 
						glTexCoord2f(1,1); 
			glVertex2f(1,1); 
			glTexCoord2f(0,1); 
			glVertex2f(-1,1); 
		}
		glEnd();
		
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		glMatrixMode(GL_MODELVIEW);
		
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D,  0);
	}
};

