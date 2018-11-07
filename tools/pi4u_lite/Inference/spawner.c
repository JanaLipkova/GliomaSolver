//#define _XOPEN_SOURCE 500
//#define _BSD_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <ftw.h>
#include <unistd.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdlib.h>

/*** HELPER ROUTINES (TBD: SOURCE) ***/
#define MY_GETTIME	1
double my_gettime()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return (double)t.tv_sec + (double)t.tv_usec*1.0E-6;
}

void parse(char *line, char **argv)
{
	while (*line != '\0') {		/* if not the end of line ....... */ 
		while (*line == ' ' || *line == '\t' || *line == '\n')
			*line++ = '\0';		/* replace white spaces with 0 */
		*argv++ = line;		/* save the argument position */
		while (*line != '\0' && *line != ' ' && 
			*line != '\t' && *line != '\n') 
			line++;	/* skip the argument until ...*/
	}
	*argv = '\0';	/* mark the end of argument list */
}


int execute_cmd(int me, char *largv[], char *dir)
{
	int rf, res = 0;

	rf = fork();
	if (rf < 0) {
		printf("spanwer(%d): fork failed!!!!\n", me); fflush(0);
		return 1;
	}

	if (rf == 0) {
		if (dir != NULL)
				chdir(dir);	/* move to the specified directory */

		//int res;
		//res =
		execvp(*largv, largv);
	}
	waitpid(rf, NULL, 0);
	return res;
}


int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
	int rv = remove(fpath);

	if (rv)
		perror(fpath);

	return rv;
}

int rmrf(char *path)
{
	return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}


int cp(const char *from, const char *to)
{
	int fd_to, fd_from;
	char buf[4096];
	ssize_t nread;
	int saved_errno;
	struct stat sb;

	fd_from = open(from, O_RDONLY);
	if (fd_from < 0)
		return -1;

	fstat(fd_from, &sb);
	if (S_ISDIR(sb.st_mode)) {	/* more supported than DT_REG */
		//printf("It is a directory!\n");
		fd_to = -1;
		goto out_error;
	}

	fd_to = open(to, O_WRONLY | O_CREAT | O_EXCL, sb.st_mode);
	if (fd_to < 0)
		goto out_error;

	while (nread = read(fd_from, buf, sizeof buf), nread > 0) {
		char *out_ptr = buf;
		ssize_t nwritten;

		do {
			nwritten = write(fd_to, out_ptr, nread);
			if (nwritten >= 0) {
				nread -= nwritten;
				out_ptr += nwritten;
			}
			else if (errno != EINTR) {
				goto out_error;
			}
		} while (nread > 0);
	}

	if (nread == 0) {
		if (close(fd_to) < 0) {
			fd_to = -1;
			goto out_error;
		}
		else {	// peh: added due to some issues on monte rosa
			fsync(fd_to);
		}
		close(fd_from);

		/* Success! */
		return 0;
	}

out_error:
	saved_errno = errno;

	close(fd_from);
	if (fd_to >= 0)
		close(fd_to);

	errno = saved_errno;
	return -1;
}

int copy_from_dir(char *name)
{
	DIR *dir;
	struct dirent *ent;
	//struct stat sb;

	dir = opendir (name);
	if (dir != NULL) {
		/* print all the files and directories within directory */
		while ((ent = readdir (dir)) != NULL) {
			//if (ent->d_type == DT_REG) {
				//printf ("%s (%d)\n", ent->d_name, ent->d_type);
				char source[256], dest[256];

				sprintf(source, "%s/%s", name, ent->d_name);
				sprintf(dest, "./%s", ent->d_name);
				cp(source, dest);
			//}

		}
		closedir (dir);
	} else {
		/* could not open directory */
		perror ("oops!");
		return 1;
	}
	return 0;
}
