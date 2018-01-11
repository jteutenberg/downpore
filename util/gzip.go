package util

import (
	"compress/gzip"
	"io"
)

//SeekableGZip is a Reader (wrapping a gzip.Reader) and adds a Seek function
//This can therefore be used as a ReadSeekCloser for sequence.SequenceSet
type SeekableGZipReader struct {
	reader *gzip.Reader
	offset int64
	buf []byte
}

func NewSeekableGZipReader(r io.Reader) *SeekableGZipReader {
	gr, _ := gzip.NewReader(r)
	return &SeekableGZipReader{reader:gr,offset:0,buf:make([]byte,65536,65536)}
}

func (s *SeekableGZipReader) Close() error {
	return s.reader.Close()
}
func (s *SeekableGZipReader) Multistream(ok bool) {
	s.reader.Multistream(ok)
}
func (s *SeekableGZipReader) Read(p []byte) (int,error) {
	n, err := s.reader.Read(p)
	s.offset += int64(n)
	return n, err
}

func (s *SeekableGZipReader) Reset(r io.Reader) error {
	s.offset = 0
	return s.reader.Reset(r)
}

func (s *SeekableGZipReader) Seek(offset int64, whence int) (int64,error) {
	if whence == io.SeekStart {
		offset -= s.offset //becomes relative to current position
	}
	length := int64(len(s.buf))
	for offset > length {
		if n, err := s.reader.Read(s.buf); err == nil || n > 0 {
			in := int64(n)
			s.offset += in
			offset -= in
		} else {
			return s.offset, err
		}
	}
	for offset > 0 {
		if n, err := s.reader.Read(s.buf[:offset]); err == nil || n > 0 {
			in := int64(n)
			s.offset += in
			offset -= in
		} else {
			return s.offset, err
		}
	}
	return s.offset, nil
}
