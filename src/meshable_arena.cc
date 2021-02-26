// -*- mode: c++; c-basic-offset: 2; indent-tabs-mode: nil -*-
// Copyright 2019 The Mesh Authors. All rights reserved.
// Use of this source code is governed by the Apache License,
// Version 2.0, that can be found in the LICENSE file.

#ifdef __linux__
#define USE_MEMFD 1
#include <linux/fs.h>
#endif
// #undef USE_MEMFD

#ifdef USE_MEMFD
#include <sys/syscall.h>
#include <unistd.h>

//#include <sys/memfd.h>
//#include <asm/unistd_64.h>
#include <sys/syscall.h>
#include <linux/memfd.h>
#endif

#include <sys/ioctl.h>

#include <algorithm>

#include "meshable_arena.h"
#include "mini_heap.h"
#include "runtime.h"

namespace mesh {

static void *arenaInstance;

static const char *const TMP_DIRS[] = {
    "/dev/shm",
    "/tmp",
};

MeshableArena::MeshableArena() : SuperHeap(), _fastPrng(internal::seed(), internal::seed()) {
  d_assert(arenaInstance == nullptr);
  arenaInstance = this;

  int fd = -1;
  if (kMeshingEnabled) {
    fd = openSpanFile(kArenaSize);
    if (fd < 0) {
      debug("mesh: opening arena file failed.\n");
      abort();
    }
  }
  _fd = fd;
  _arenaBegin = SuperHeap::map(kArenaSize, kMapShared, fd);
  _mhIndex = reinterpret_cast<atomic<MiniHeapID> *>(SuperHeap::malloc(indexSize()));

  hard_assert(_arenaBegin != nullptr);
  hard_assert(_mhIndex != nullptr);

  if (kAdviseDump) {
    madvise(_arenaBegin, kArenaSize, MADV_DONTDUMP);
  }

  // debug("MeshableArena(%p): fd:%4d\t%p-%p\n", this, fd, _arenaBegin, arenaEnd());

  // TODO: move this to runtime
  atexit(staticAtExit);
  pthread_atfork(staticPrepareForFork, staticAfterForkParent, staticAfterForkChild);
}

char *uintToStr(char *dst, uint32_t i) {
  constexpr size_t maxLen = sizeof("4294967295") + 1;
  char buf[maxLen];
  memset(buf, 0, sizeof(buf));

  char *digit = buf + maxLen - 2;
  // capture the case where i == 0
  *digit = '0';
  while (i > 0) {
    hard_assert(reinterpret_cast<uintptr_t>(digit) >= reinterpret_cast<uintptr_t>(buf));
    const char c = '0' + (i % 10);
    *digit = c;
    digit--;
    i /= 10;
  }
  if (*digit == '\0') {
    digit++;
  }

  return strcat(dst, digit);
}

char *MeshableArena::openSpanDir(int pid) {
  constexpr size_t buf_len = 128;

  for (auto tmpDir : TMP_DIRS) {
    for (size_t i = 0; i < 1024; i++) {
      char buf[buf_len];
      memset(buf, 0, buf_len);

      // on some platforms snprintf actually calls out to malloc,
      // despite us passing in a reasonable buffer.  Since what we're doing is
      // reasonably simple, just build the path ourselves to avoid this.
      char *next = buf;
      hard_assert(strlen(tmpDir) < buf_len);
      next = strcat(next, tmpDir);
      next = strcat(next, "/alloc-mesh-");
      next = uintToStr(next, pid);
      next = strcat(next, ".");
      next = uintToStr(next, i);

      // ensure we haven't overflown our buffer
      hard_assert(reinterpret_cast<uintptr_t>(next) <= reinterpret_cast<uintptr_t>(buf) + buf_len);

      int result = mkdir(buf, 0755);
      if (result != 0) {
        if (errno == EEXIST) {
          // we will get EEXIST if we have re-execed -- we need to use a
          // new directory because we could have dropped privileges in
          // the meantime.
          continue;
        } else {
          // otherwise it is likely that the parent tmp directory
          // doesn't exist or we don't have permissions in it.
          break;
        }
      }

      char *spanDir = reinterpret_cast<char *>(internal::Heap().malloc(strlen(buf) + 1));
      strcpy(spanDir, buf);
      return spanDir;
    }
  }

  return nullptr;
}

void MeshableArena::expandArena(size_t minPagesAdded) {
  // debug("expandArena : %d, end=%d\n", minPagesAdded, _end);
  const size_t pageCount = std::max(minPagesAdded, kMinArenaExpansion);

  Span expansion(_end, pageCount);
  _end += pageCount;

  if (unlikely(_end >= kArenaSize / kPageSize)) {
    debug("Mesh: arena exhausted: current arena size is %.1f GB; recompile with larger arena size.",
          kArenaSize / 1024.0 / 1024.0 / 1024.0);
    abort();
  }

  // for(size_t i = 0; i < kSpanClassCount; ++i) {
  //   if(!_dirty[i].empty()) {
  //    debug("_dirty class[%d] = %d\n", i, _dirty[i].size());
  //   }
  // }

  // for(size_t i = 0; i < kSpanClassCount; ++i) {
  //   if(!_clean[i].empty()) {
  //    debug("_clean class[%d] = %d\n", i, _clean[i].size());
  //   }
  // }

  // for(size_t i = 0; i < kSpanClassCount; ++i) {
  //   if(!_dirty[i].empty()) {
  //     debug("_dirty[%d]  size = %d", i, _dirty[i].size());
  //   }
  // }

  // for(size_t i = 0; i < kSpanClassCount; ++i) {
  //   if(!_clean[i].empty()) {
  //     debug("_clean[%d]  size = %d", i, _clean[i].size());
  //   }
  // }

  // for(auto& s : _clean[kSpanClassCount-1]) {
  //   debug("clean(%d, %d)", s.offset, s.length);
  // }
  // debug("expandArena (minPagesAdded=%d) %d, %d,  dirty : %d\n", minPagesAdded, expansion.offset, expansion.length,
  // _dirtyPageCount);

  if (_isCOWRunning) {
    trackCOWed(expansion);
    resetSpanMapping(expansion);
  }
  _clean[expansion.spanClass()].push_back(expansion);
  // debug("expandArena : %d, end=%d\n", minPagesAdded, _end);
}

bool MeshableArena::findPagesInner(internal::vector<Span> freeSpans[kSpanClassCount], const size_t i,
                                   const size_t pageCount, Span &result) {
  internal::vector<Span> &spanList = freeSpans[i];
  if (spanList.empty())
    return false;

  size_t oldLen = spanList.size();

  if (i == kSpanClassCount - 1) {
    // the final span class contains (and is the only class to
    // contain) variable-size spans, so we need to make sure we
    // search through all candidates in this case.
    for (int64_t j = spanList.size() - 1; j > 0; --j) {
      if (spanList[j].length >= pageCount) {
        std::swap(spanList[j], spanList.back());
        break;
      }
    }

    // check that we found something in the above loop. this would be
    // our last loop iteration anyway
    if (spanList.back().length < pageCount) {
      return false;
    }
  }

  Span span = spanList.back();
  spanList.pop_back();

#ifndef NDEBUG
  d_assert_msg(oldLen == spanList.size() + 1, "pageCount:%zu,%zu -- %zu/%zu", pageCount, i, oldLen, spanList.size());
  for (size_t j = 0; j < spanList.size(); j++) {
    d_assert(spanList[j] != span);
  }
#endif

  // this invariant should be maintained
  d_assert(span.length >= i + 1);
  d_assert(span.length >= pageCount);

  // put the part we don't need back in the reuse pile
  Span rest = span.splitAfter(pageCount);

  auto ci = rest.spanClass();
  if (!rest.empty()) {
    freeSpans[ci].push_back(rest);
    moveBiggerTofirst(freeSpans[ci]);
  }
  d_assert(span.length == pageCount);

  result = span;
  return true;
}

bool MeshableArena::findPagesInnerFast(internal::vector<Span> freeSpans[kSpanClassCount], const size_t i,
                                       const size_t pageCount, Span &result) {
  internal::vector<Span> &spanList = freeSpans[i];
  if (unlikely(spanList.empty()))
    return false;

  Span span = spanList.back();
  spanList.pop_back();

  // this invariant should be maintained
  d_assert(span.length >= i + 1);
  d_assert(span.length >= pageCount);

  // put the part we don't need back in the reuse pile
  Span rest = span.splitAfter(pageCount);
  if (!rest.empty()) {
    freeSpans[rest.spanClass()].push_back(rest);
    moveBiggerTofirst(freeSpans[rest.spanClass()]);
  }
  d_assert(span.length == pageCount);

  result = span;
  return true;
}

bool MeshableArena::findPages(const size_t pageCount, Span &result, internal::PageType &type) {
  auto targetClass = Span(0, pageCount).spanClass();

  // search the fix length class first
  for (size_t i = targetClass; i < kSpanClassCount - 1; ++i) {
    if (findPagesInnerFast(_dirty, i, pageCount, result)) {
      type = internal::PageType::Dirty;
      _dirtyPageCount -= result.length;
      return true;
    }

    if (findPagesInnerFast(_clean, i, pageCount, result)) {
      type = internal::PageType::Clean;
      return true;
    }
  }

  // Search through all dirty spans first.  We don't worry about
  // fragmenting dirty pages, as being able to reuse dirty pages means
  // we don't increase RSS.
  if (findPagesInner(_dirty, kSpanClassCount - 1, pageCount, result)) {
    type = internal::PageType::Dirty;
    _dirtyPageCount -= result.length;
    return true;
  }

  // if no dirty pages are available, search clean pages.  An allocated
  // clean page (once it is written to) means an increased RSS.
  if (findPagesInner(_clean, kSpanClassCount - 1, pageCount, result)) {
    type = internal::PageType::Clean;
    return true;
  }

  return false;
}

Span MeshableArena::reservePages(const size_t pageCount, const size_t pageAlignment) {
  d_assert(pageCount >= 1);

  internal::PageType flags(internal::PageType::Unknown);
  Span result(0, 0);
  auto ok = findPages(pageCount, result, flags);
  if (!ok) {
    tryAndSendToFree(new internal::FreeCmd(internal::FreeCmd::FLUSH));
    getSpansFromBg(true);
    ok = findPages(pageCount, result, flags);
    if (!ok) {
      expandArena(pageCount);
      ok = findPages(pageCount, result, flags);
      hard_assert(ok);
    }
  }

  d_assert(!result.empty());
  d_assert(flags != internal::PageType::Unknown);

  if (unlikely(pageAlignment > 1 && ((ptrvalFromOffset(result.offset) / kPageSize) % pageAlignment != 0))) {
    freeSpan(result, flags);
    // recurse once, asking for enough extra space that we are sure to
    // be able to find an aligned offset of pageCount pages within.
    result = reservePages(pageCount + 2 * pageAlignment, 1);

    const size_t alignment = pageAlignment * kPageSize;
    const uintptr_t alignedPtr = (ptrvalFromOffset(result.offset) + alignment - 1) & ~(alignment - 1);
    const auto alignedOff = offsetFor(reinterpret_cast<void *>(alignedPtr));
    d_assert(alignedOff >= result.offset);
    d_assert(alignedOff < result.offset + result.length);
    const auto unwantedPageCount = alignedOff - result.offset;
    auto alignedResult = result.splitAfter(unwantedPageCount);
    d_assert(alignedResult.offset == alignedOff);
    freeSpan(result, flags);
    const auto excess = alignedResult.splitAfter(pageCount);
    freeSpan(excess, flags);
    result = alignedResult;
  }

  return result;
}

template <typename Func>
static void forEachFree(const internal::vector<Span> freeSpans[kSpanClassCount], const Func func) {
  for (size_t i = 0; i < kSpanClassCount; i++) {
    if (freeSpans[i].empty())
      continue;

    for (size_t j = 0; j < freeSpans[i].size(); j++) {
      auto span = freeSpans[i][j];
      func(span);
    }
  }
}

internal::RelaxedBitmap MeshableArena::allocatedBitmap(bool includeDirty) const {
  internal::RelaxedBitmap bitmap(_end);

  // we can build up a bitmap of in-use pages here by looking at the
  // arena start and end addresses (to compute the number of
  // bits/pages), set all bits to 1, then iterate through our _clean
  // and _dirty lists unsetting pages that aren't in use.

  bitmap.setAll(_end);

  auto unmarkPages = [&](const Span &span) {
    for (size_t k = 0; k < span.length; k++) {
#ifndef NDEBUG
      if (!bitmap.isSet(span.offset + k)) {
        debug("arena: bit %zu already unset 1 (%zu/%zu)\n", k, span.offset, span.length);
      }
#endif
      bitmap.unset(span.offset + k);
    }
  };

  if (includeDirty)
    forEachFree(_dirty, unmarkPages);
  forEachFree(_clean, unmarkPages);

  return bitmap;
}

char *MeshableArena::pageAlloc(Span &result, size_t pageCount, size_t pageAlignment) {
  if (pageCount == 0) {
    return nullptr;
  }

  d_assert(_arenaBegin != nullptr);

  d_assert(pageCount >= 1);
  d_assert(pageCount < std::numeric_limits<Length>::max());

  auto span = reservePages(pageCount, pageAlignment);
  d_assert(isAligned(span, pageAlignment));

  d_assert(contains(ptrFromOffset(span.offset)));
#ifndef NDEBUG
  if (_mhIndex[span.offset].load().hasValue()) {
    mesh::debug("----\n");
    auto mh = reinterpret_cast<MiniHeap *>(miniheapForArenaOffset(span.offset));
    mh->dumpDebug();
  }
#endif

  char *ptr = reinterpret_cast<char *>(ptrFromOffset(span.offset));

  result = span;
  return ptr;
}

void MeshableArena::free(void *ptr, size_t sz, internal::PageType type) {
  if (unlikely(!contains(ptr))) {
    debug("invalid free of %p/%zu", ptr, sz);
    return;
  }
  d_assert(sz > 0);

  d_assert(sz / kPageSize > 0);
  d_assert(sz % kPageSize == 0);

  const Span span(offsetFor(ptr), sz / kPageSize);
  freeSpan(span, type);
}

static size_t flushSpansByOffset(internal::vector<Span> freeSpans[kSpanClassCount], internal::vector<Span> &flushSpans,
                                 Offset offsetBegin, size_t needFreeCount) {
  size_t freeCount = 0;
  Offset offsetEnd = offsetBegin + needFreeCount;

  flushSpans.reserve(needFreeCount);

  for (size_t i = 0; i < kSpanClassCount; ++i) {
    auto &spans = freeSpans[i];

    if (spans.empty())
      continue;

    internal::vector<Span> rest;
    rest.reserve(spans.size());

    for (auto &span : spans) {
      if (offsetBegin <= span.offset && span.offset < offsetEnd) {
        flushSpans.emplace_back(span);
        freeCount += span.length;
      } else {
        rest.emplace_back(span);
      }
    }
    spans.swap(rest);
  }

  return freeCount;
}

struct {
  bool operator()(const Span &a, const Span &b) const {
    return a.length > b.length;
  }
} customLess;

void MeshableArena::getSpansFromBg(bool wait) {
  bool needSort = false;
  bool gotOne = false;

  constexpr unsigned int spin_loops = 16u, spins = 16u;
  unsigned int loop_count = 0;

  while (true) {
    size_t pageCount = 0;
    internal::FreeCmd *preCommand = runtime().getReturnCmdFromBg();
    ++loop_count;

    if (preCommand) {
      // debug("getSpansFromBg = %d\n", preDirtySpans->size());
      // all add to the mark spans
      if (preCommand->cmd == internal::FreeCmd::FLUSH) {
        for (auto &s : preCommand->spans) {
#ifndef NDEBUG
          if (_isCOWRunning) {
            for (auto i = s.offset; i < s.length; ++i) {
              d_assert(_cowBitmap.isSet(s.offset));
            }
          }
#endif
          _clean[s.spanClass()].emplace_back(s);
          pageCount += s.length;
        }

        if (preCommand->spans.size() > 0) {
          needSort = true;
        }
        gotOne = true;
      } else {
        hard_assert(false);
      }
      // debug("getSpansFromBg got %d spans -  %d page from backgroud.\n", preCommand->spans.size(), pageCount);
      delete preCommand;
    } else {
      if (wait && !gotOne && runtime().freeThreadRunning()) {
        if (loop_count < spin_loops) {
          for (unsigned int j = 0; j < spins; ++j) {
            cpupause();
          }
        } else {
          std::this_thread::yield();
        }
        continue;
      } else {
        break;
      }
    }
  }

  if (needSort) {
    std::sort(_clean[kSpanClassCount - 1].begin(), _clean[kSpanClassCount - 1].end(), customLess);
  }
  // debug("getSpansFromBg after sort last");
  // for(size_t i = 0; i < kSpanClassCount; ++i) {
  //   if(!_dirty[i].empty()) {
  //     debug("_dirty[%d]  size = %d", i, _dirty[i].size());
  //   }
  // }

  // for(size_t i = 0; i < kSpanClassCount; ++i) {
  //   if(!_clean[i].empty()) {
  //     debug("_clean[%d]  size = %d", i, _clean[i].size());
  //   }
  // }

  // for(auto& s : _clean[kSpanClassCount-1]) {
  //   debug("clean(%d, %d)", s.offset, s.length);
  // }
  // debug("getSpansFromBg end");
}

void MeshableArena::tryAndSendToFree(internal::FreeCmd *fCommand) {
  auto &rt = runtime();
  bool ok = rt.sendFreeCmd(fCommand);
  while (!ok) {
    getSpansFromBg();
    ok = rt.sendFreeCmd(fCommand);
  }
}

void MeshableArena::dumpSpans() {
  size_t dirty = 0;
  size_t clean = 0;
  for (size_t i = 0; i < kSpanClassCount; ++i) {
    if (_dirty[i].size() || _clean[i].size()) {
      debug("MeshInfo spanClass:%-3lu  dirty:%6zu, clean:%6zu", i, _dirty[i].size(), _clean[i].size());
      dirty += _dirty[i].size();
      clean += _clean[i].size();
    }
  }
  debug("MeshInfo span summary reset: %zu, dirty: %zu, clean: %zu", _toReset.size(), dirty, clean);
}

void MeshableArena::partialScavenge() {
  size_t needFreeCount = 0;

  if (_dirtyPageCount > kMaxDirtyPageThreshold) {
    needFreeCount = _end / 5;
  }

  internal::FreeCmd *freeCommand = new internal::FreeCmd(internal::FreeCmd::FREE_DIRTY_PAGE);
  size_t freeCount = flushSpansByOffset(_dirty, freeCommand->spans, _lastFlushBegin, needFreeCount);

  _dirtyPageCount -= freeCount;

  _lastFlushBegin += needFreeCount;
  if (_lastFlushBegin >= _end) {
    _lastFlushBegin = 0;
  }

  tryAndSendToFree(freeCommand);
}

void MeshableArena::scavenge(bool force) {
  if (!force && _dirtyPageCount < kMinDirtyPageThreshold) {
    return;
  }

  internal::FreeCmd *unmapCommand = new internal::FreeCmd(internal::FreeCmd::UNMAP_PAGE);

  auto markPages = [&](const Span &span) {
    // debug("arena:  (%zu/%zu) \n", span.offset, span.length);
    unmapCommand->spans.emplace_back(span);
  };

  // first, untrack the spans in the meshed bitmap and mark them in
  // the (method-local) unallocated bitmap
  std::for_each(_toReset.begin(), _toReset.end(), [&](const Span &span) {
    markPages(span);
    // resetSpanMapping(span);
  });

  // now that we've finally reset to identity all delayed-reset
  // mappings, empty the list
  // debug("_toReset size: %d", _toReset.size());
  _toReset.clear();

  tryAndSendToFree(unmapCommand);

  internal::FreeCmd *freeCommand = new internal::FreeCmd(internal::FreeCmd::FREE_DIRTY_PAGE);

  size_t needFreeCount = _end / 5;

  size_t freeCount = flushSpansByOffset(_dirty, freeCommand->spans, _lastFlushBegin, needFreeCount);

  _dirtyPageCount -= freeCount;

  tryAndSendToFree(freeCommand);

  internal::FreeCmd *cleanCommand = new internal::FreeCmd(internal::FreeCmd::CLEAN_PAGE);
  freeCount = flushSpansByOffset(_clean, cleanCommand->spans, _lastFlushBegin, needFreeCount);

  _lastFlushBegin += needFreeCount;
  if (_lastFlushBegin >= _end) {
    _lastFlushBegin = 0;
  }

  tryAndSendToFree(cleanCommand);
  // debug("FreeCmd::CLEAN_PAGE");

  tryAndSendToFree(new internal::FreeCmd(internal::FreeCmd::FLUSH));

  getSpansFromBg();
}

void MeshableArena::freePhys(const Span &span) {
  auto ptr = ptrFromOffset(span.offset);
  auto sz = span.byteLength();
  freePhys(ptr, sz);
}

void MeshableArena::freePhys(void *ptr, size_t sz) {
  d_assert(contains(ptr));
  d_assert(sz > 0);

  d_assert(sz / CPUInfo::PageSize > 0);
  d_assert(sz % CPUInfo::PageSize == 0);

  // we madvise(MADV_DONTNEED) elsewhere; this function is only needed
  // when our heap is a shared mapping
  if (!kMeshingEnabled) {
    return;
  }

  const off_t off = reinterpret_cast<char *>(ptr) - reinterpret_cast<char *>(_arenaBegin);
#ifndef __APPLE__
  int result = fallocate(_fd, FALLOC_FL_PUNCH_HOLE | FALLOC_FL_KEEP_SIZE, off, sz);
  d_assert_msg(result == 0, "result(fd %d): %d errno %d (%s)\n", _fd, result, errno, strerror(errno));
#else
#warning macOS version of fallocate goes here
  fstore_t store = {F_ALLOCATECONTIG, F_PEOFPOSMODE, 0, (long long)sz, 0};
  int result = fcntl(_fd, F_PREALLOCATE, &store);
  if (result == -1) {
    // try and allocate space with fragments
    store.fst_flags = F_ALLOCATEALL;
    result = fcntl(_fd, F_PREALLOCATE, &store);
  }
  // if (result != -1) {
  //    result = ftruncate(_fd, off+sz);
  // }
  d_assert(result == 0);
#endif
}

void MeshableArena::beginMesh(void *keep, void *remove, size_t sz) {
  int r = mprotect(remove, sz, PROT_READ);
  hard_assert(r == 0);
}

void MeshableArena::finalizeMesh(void *keep, void *remove, size_t sz) {
  // debug("keep: %p, remove: %p\n", keep, remove);
  const auto keepOff = offsetFor(keep);
  // const auto removeOff = offsetFor(remove);

  const size_t pageCount = sz / kPageSize;

  hard_assert(pageCount < std::numeric_limits<Length>::max());
  // const Span removedSpan{removeOff, static_cast<Length>(pageCount)};

  void *ptr = mmap(remove, sz, HL_MMAP_PROTECTION_MASK, kMapShared | MAP_FIXED, _fd, keepOff * kPageSize);
  hard_assert_msg(ptr != MAP_FAILED, "mesh remap failed: %d", errno);
}

int MeshableArena::openShmSpanFile(size_t sz) {
  constexpr size_t buf_len = 64;
  char buf[buf_len];
  memset(buf, 0, buf_len);

  _spanDir = openSpanDir(getpid());
  d_assert(_spanDir != nullptr);

  char *next = strcat(buf, _spanDir);
  strcat(next, "/XXXXXX");

  int fd = mkstemp(buf);
  if (fd < 0) {
    debug("mkstemp: %d (%s)\n", errno, strerror(errno));
    abort();
  }

  // we only need the file descriptors, not the path to the file in the FS
  int err = unlink(buf);
  if (err != 0) {
    debug("unlink: %d\n", errno);
    abort();
  }

  // TODO: see if fallocate makes any difference in performance
  err = ftruncate(fd, sz);
  if (err != 0) {
    debug("ftruncate: %d\n", errno);
    abort();
  }

  // if a new process gets exec'ed, ensure our heap is completely freed.
  err = fcntl(fd, F_SETFD, FD_CLOEXEC);
  if (err != 0) {
    debug("fcntl: %d\n", errno);
    abort();
  }

  return fd;
}

#ifdef USE_MEMFD
static int sys_memfd_create(const char *name, unsigned int flags) {
  return syscall(__NR_memfd_create, name, flags);
}

int MeshableArena::openSpanFile(size_t sz) {
  errno = 0;
  int fd = sys_memfd_create("mesh_arena", MFD_CLOEXEC);
  // the call to memfd failed -- fall back to opening a shm file
  if (fd < 0) {
    return openShmSpanFile(sz);
  }

  int err = ftruncate(fd, sz);
  if (err != 0) {
    debug("ftruncate: %d\n", errno);
    abort();
  }

  return fd;
}
#else
int MeshableArena::openSpanFile(size_t sz) {
  return openShmSpanFile(sz);
}
#endif  // USE_MEMFD

void MeshableArena::staticAtExit() {
  d_assert(arenaInstance != nullptr);
  if (arenaInstance != nullptr)
    reinterpret_cast<MeshableArena *>(arenaInstance)->exit();
}

void MeshableArena::staticPrepareForFork() {
  d_assert(arenaInstance != nullptr);
  reinterpret_cast<MeshableArena *>(arenaInstance)->prepareForFork();
}

void MeshableArena::staticAfterForkParent() {
  d_assert(arenaInstance != nullptr);
  reinterpret_cast<MeshableArena *>(arenaInstance)->afterForkParent();
}

void MeshableArena::staticAfterForkChild() {
  d_assert(arenaInstance != nullptr);
  reinterpret_cast<MeshableArena *>(arenaInstance)->afterForkChild();
}

void MeshableArena::prepareForFork() {
  if (!kMeshingEnabled) {
    return;
  }

  runtime().heap().lock();
  runtime().lock();
  internal::Heap().lock();

  // block here until the COW is finished.
  while (_isCOWRunning) {
    moveRemainPages();
  }

  tryAndSendToFree(new internal::FreeCmd(internal::FreeCmd::FLUSH));
  getSpansFromBg(true);

  internal::Heap().lock();

  _isCOWRunning = true;
  hard_assert(_lastCOW == 0);

  _preSpanDir = _spanDir;
  _prefd = _fd;

  _fd = openSpanFile(kArenaSize);

  struct stat fileinfo;
  memset(&fileinfo, 0, sizeof(fileinfo));
  fstat(_fd, &fileinfo);
  d_assert(fileinfo.st_size >= 0 && (size_t)fileinfo.st_size == kArenaSize);

  // debug("%d: prepare fork", getpid());
  int r = mprotect(_arenaBegin, kArenaSize, PROT_READ);
  hard_assert(r == 0);
}

void MeshableArena::afterForkParentAndChild() {
  // remap the large area, from end to the big end
  size_t hasnt_use = kArenaSize - _end * kPageSize;
  void *address = (void *)((size_t)_arenaBegin + _end * kPageSize);
  size_t address_size = hasnt_use;
  size_t address_offset = _end * kPageSize;
  _COWend = _end;

<<<<<<< HEAD
  void *ptr = mmap(address, address_size, HL_MMAP_PROTECTION_MASK, kMapShared | MAP_FIXED, _fd, address_offset);
  hard_assert_msg(ptr == address, "map failed: %d, addr=%p, %u, %u", errno, address, address_size, address_offset);
#ifndef NDEBUG
  debug("afterForkParentAndChild remap %d: errno=%d, addr=%p, %u, %u", getpid(), errno, address, address_size,
        address_offset);
#endif
=======
  internal::Heap().unlock();

  close(_forkPipe[1]);
>>>>>>> f66adeb9cfd00c463e055ce6b33b64625ed26f5b

  if (kAdviseDump) {
    madvise(address, address_size, MADV_DONTDUMP);
  }

  tryAndSendToFree(new internal::FreeCmd(internal::FreeCmd::FLUSH));
  getSpansFromBg(true);

  // remap all the clean spans
  size_t count = 0;
  for (size_t i = 0; i < kSpanClassCount; ++i) {
    for (auto &span : _clean[i]) {
      trackCOWed(span);
      resetSpanMapping(span);
      ++count;
    }
  }

#ifndef NDEBUG
  debug("afterForkParentAndChild %d: resetSpanMapping %u clean spans", getpid(), count);
#endif
  count = 0;
  for (size_t i = 0; i < kSpanClassCount; ++i) {
    for (auto &span : _dirty[i]) {
      trackCOWed(span);
      resetSpanMapping(span);
    }
    ++count;
  }
#ifndef NDEBUG
  debug("afterForkParentAndChild %d: resetSpanMapping %u dirty spans", getpid(), count);
#endif
  count = 0;
  for (auto &span : _toReset) {
    trackCOWed(span);
    resetSpanMapping(span);
    ++count;
  }
#ifndef NDEBUG
  debug("afterForkParentAndChild %d: resetSpanMapping %u _toReset spans", getpid(), count);
#endif
}

void MeshableArena::afterForkParent() {
  internal::Heap().unlock();

  if (kMeshingEnabled) {
    afterForkParentAndChild();
  }

  runtime().unlock();
  runtime().heap().unlock();
}

void MeshableArena::doAfterForkChild() {
  afterForkChild();
}

void MeshableArena::afterForkChild() {
  if (runtime().pid() == getpid()) {
    return;
  }

  runtime().updatePid();
  debug("afterForkChild pid=%u", getpid());

  runtime().setFreeThreadRunning(false);
  _needCOWScan = false;

  if (!kMeshingEnabled) {
    return;
  }

  close(_fd);

  _fd = openSpanFile(kArenaSize);

  struct stat fileinfo;
  memset(&fileinfo, 0, sizeof(fileinfo));
  fstat(_fd, &fileinfo);
  d_assert(fileinfo.st_size >= 0 && (size_t)fileinfo.st_size == kArenaSize);

  internal::Heap().unlock();

  afterForkParentAndChild();

  runtime().unlock();
  runtime().heap().unlock();
}

bool MeshableArena::moveMiniHeapToNewFile(MiniHeap *mh, void *ptr) {
  // debug("moveMiniHeapToNewFile %d: mh=%p  addr=%p\n", getpid(), mh, ptr);

  MiniHeap *leader_mh = mh->meshedLeader();

  if (leader_mh != nullptr) {
    debug("moveMiniHeapToNewFile %d: mh=%p  leader=%p\n", getpid(), mh, leader_mh);
  } else {
    leader_mh = mh;
  }

  const auto sz = leader_mh->spanSize();
  auto keep = reinterpret_cast<void *>(leader_mh->getSpanStart(arenaBegin()));
  const auto keepOff = offsetFor(keep);

  const auto &span = leader_mh->span();

  // only copy the phys page once
  for (size_t i = 0; i < span.length; ++i) {
    if (_cowBitmap.isSet(span.offset + i)) {
      hard_assert_msg(i == 0, "moveMiniHeapToNewFile %d: already COW span offset=%u, length=%u", getpid(),
                      span.offset + i, span.length);
      debug("moveMiniHeapToNewFile %d: trigger doulbe move, just return and try again !!!", getpid(), span.offset + i,
            span.length);
      return true;
    }
    d_assert(span.offset == keepOff);
  }

  // copy the phys pages to the new file
  auto copyoff = keepOff * kPageSize;
  auto copysize = span.length * kPageSize;
  auto result = internal::copyFile(_fd, _prefd, copyoff, copysize);

  hard_assert_msg((size_t)result == copysize, "internal::copyFile err rc=%d, %u, %u, %u, %u\n", result, _fd, _prefd,
                  copyoff, copysize);

  // redo all the mapping
  const auto base = leader_mh;

  base->forEachMeshed([&](MiniHeap *mh) {
    trackCOWed(mh->span());

    keep = reinterpret_cast<void *>(mh->getSpanStart(arenaBegin()));

    void *remap_ptr = mmap(keep, sz, HL_MMAP_PROTECTION_MASK, kMapShared | MAP_FIXED, _fd, keepOff * kPageSize);

    hard_assert_msg(remap_ptr != MAP_FAILED, "mesh remap failed: %d", errno);
    // debug("moveMiniHeapToNewFile %d:  remap addr=%p, size=%d, mh=%p\n", getpid(), remap_ptr, sz, mh);
    return false;
  });

  return true;
}

void MeshableArena::moveRemainPages() {
#ifndef NDEBUG
  debug("moveRemainPages checking : %u / %u, %u", _lastCOW, _COWend, _end);
#endif

  auto check_begin = _lastCOW;

  size_t off;
  for (off = check_begin; off < check_begin + kMaxCOWPage && off < _COWend;) {
    d_assert(off != _COWend);
    MiniHeap *mh = reinterpret_cast<MiniHeap *>(miniheapForArenaOffset(off));

    if (_cowBitmap.isSet(off)) {
      if (mh) {
#ifndef NDEBUG
        for (auto i = off; i < mh->span().length; ++i) {
          d_assert(_cowBitmap.isSet(off));
        }
#endif
        MiniHeap *leader_mh = mh->meshedLeader();
        d_assert_msg(mh->span().offset <= off && off < mh->span().offset + mh->span().length,
                     "mh->span(%u, %u)  off=%u  leader=%p", mh->span().offset, mh->span().length, off, leader_mh);
        off += mh->span().length;
      } else {
        ++off;
      }
      continue;
    } else {
      if (mh) {
        moveMiniHeapToNewFile(mh, nullptr);
        d_assert(_cowBitmap.isSet(off));
        d_assert(off == mh->span().offset || off == mh->meshedLeader()->span().offset);
        off += mh->span().length;
      } else {
        tryAndSendToFree(new internal::FreeCmd(internal::FreeCmd::FLUSH));
        getSpansFromBg(true);

        bool found = false;
        Offset slen = 0;
        for (size_t i = 0; i < kSpanClassCount; ++i) {
          for (auto &span : _clean[i]) {
            if (span.offset == off) {
              found = true;
              slen = span.length;
              trackCOWed(span);
              resetSpanMapping(span);
#ifndef NDEBUG
              debug("found in _clean len=%u", slen);
#endif
            }
          }
        }

        for (size_t i = 0; i < kSpanClassCount; ++i) {
          for (auto &span : _dirty[i]) {
            if (span.offset == off) {
              found = true;
              slen = span.length;
              trackCOWed(span);
              resetSpanMapping(span);
#ifndef NDEBUG
              debug("found in _dirty len=%u", slen);
#endif
            }
          }
        }

        for (auto &span : _toReset) {
          if (span.offset == off) {
            found = true;
            slen = span.length;
            trackCOWed(span);
            resetSpanMapping(span);
#ifndef NDEBUG
            debug("found in _toReset len=%u", slen);
#endif
          }
        }
        d_assert(found);

        off += slen;
      }
      continue;
    }
  }

  _lastCOW = off;

  if (off >= _COWend) {
#ifndef NDEBUG
    debug("moveRemainPages finished < %u, %u, %u", off, _COWend, _end);
    for (size_t j = 0; j < _COWend; ++j) {
      if (!_cowBitmap.isSet(j)) {
        tryAndSendToFree(new internal::FreeCmd(internal::FreeCmd::FLUSH));
        getSpansFromBg(true);

        bool found = false;
        Offset slen = 0;
        for (size_t i = 0; i < kSpanClassCount; ++i) {
          for (auto &span : _clean[i]) {
            if (span.offset == off) {
              found = true;
              slen = span.length;
              debug("found in _clean len=%u", slen);
            }
          }
        }

        for (size_t i = 0; i < kSpanClassCount; ++i) {
          for (auto &span : _dirty[i]) {
            if (span.offset == off) {
              found = true;
              slen = span.length;
              debug("found in _dirty len=%u", slen);
            }
          }
        }

        for (auto &span : _toReset) {
          if (span.offset == off) {
            found = true;
            slen = span.length;
            debug("found in _toReset len=%u", slen);
          }
        }
        if (!found) {
          _lastCOW = j;
          // d_assert_msg(false, "not found off= %u", j);
          return;
        }
        j += slen;
      }
    }
    debug("moveRemainPages finished %u, %u, %u >", off, _COWend, _end);
#endif

    _lastCOW = 0;
    _COWend = 0;
    _isCOWRunning = false;
    close(_prefd);
    _prefd = -1;
    _cowBitmap.clear();
  }
}

}  // namespace mesh
