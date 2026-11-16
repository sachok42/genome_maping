#include <bits/stdc++.h>
#include <zlib.h>
using namespace std;

// sa-is is a suffix array construction algorythm

vector<int> sa_is(const vector<int>& text, int alphabet_max) {
    int length = text.size();

    if (length == 0) return {};
    if (length == 1) return {0};

    if (length == 2) {
        if (text[0] < text[1]) return {0, 1};
        return {1, 0};
    }

    if (length < 40) {
        vector<int> suffix_array(length);
        iota(suffix_array.begin(), suffix_array.end(), 0);

        sort(suffix_array.begin(), suffix_array.end(),
             [&](int first, int second) {
            while (first < length && second < length) {
                if (text[first] != text[second])
                    return text[first] < text[second];

                first++;
                second++;
            }

            return first == length && second != length;
        });

        return suffix_array;
    }

    vector<int> suffix_array(length);
    vector<bool> is_s_type(length);

    for (int i = length - 2; i >= 0; i--) {
        is_s_type[i] = (text[i] == text[i + 1]
            ? is_s_type[i + 1]
            : text[i] < text[i + 1]);
    }

    vector<int> bucket_left(alphabet_max + 2, 0);
    vector<int> bucket_right(alphabet_max + 2, 0);

    for (int i = 0; i < length; i++) {
        if (!is_s_type[i])
            bucket_right[text[i]]++;
        else
            bucket_left[text[i] + 1]++;
    }

    for (int i = 0; i <= alphabet_max; i++) {
        bucket_right[i] += bucket_left[i];

        if (i < alphabet_max)
            bucket_left[i + 1] += bucket_right[i];
    }

    auto induce = [&](const vector<int>& lms_positions) {
        fill(suffix_array.begin(), suffix_array.end(), -1);

        vector<int> bucket_positions = bucket_right;

        for (int lms_position : lms_positions) {
            if (lms_position < length)
                suffix_array[bucket_positions[text[lms_position]]++] =
                    lms_position;
        }

        bucket_positions = bucket_left;
        suffix_array[bucket_positions[text[length - 1]]++] = length - 1;

        for (int i = 0; i < length; i++) {
            int suffix_position = suffix_array[i];

            if (suffix_position >= 1 &&
                !is_s_type[suffix_position - 1]) {
                int previous_position = suffix_position - 1;

                suffix_array[bucket_positions[text[previous_position]]++] =
                    previous_position;
            }
        }

        bucket_positions = bucket_left;

        for (int i = length - 1; i >= 0; i--) {
            int suffix_position = suffix_array[i];

            if (suffix_position >= 1 &&
                is_s_type[suffix_position - 1]) {
                int previous_position = suffix_position - 1;

                suffix_array[--bucket_positions[
                    text[previous_position] + 1
                ]] = previous_position;
            }
        }
    };

    vector<int> lms_map(length + 1, -1);
    int lms_count = 0;

    for (int i = 1; i < length; i++) {
        if (!is_s_type[i - 1] && is_s_type[i])
            lms_map[i] = lms_count++;
    }

    vector<int> lms_positions;
    lms_positions.reserve(lms_count);

    for (int i = 1; i < length; i++) {
        if (!is_s_type[i - 1] && is_s_type[i])
            lms_positions.push_back(i);
    }

    induce(lms_positions);

    if (lms_count > 0) {
        vector<int> sorted_lms_positions;
        sorted_lms_positions.reserve(lms_count);

        for (int suffix_position : suffix_array) {
            if (suffix_position >= 0 &&
                lms_map[suffix_position] != -1) {
                sorted_lms_positions.push_back(suffix_position);
            }
        }

        vector<int> reduced_text(lms_count);
        int reduced_alphabet_max = 0;

        reduced_text[lms_map[sorted_lms_positions[0]]] = 0;

        for (int i = 1; i < lms_count; i++) {
            int previous_lms = sorted_lms_positions[i - 1];
            int current_lms = sorted_lms_positions[i];

            int previous_end = (lms_map[previous_lms] + 1 < lms_count)
                ? lms_positions[lms_map[previous_lms] + 1]
                : length;

            int current_end = (lms_map[current_lms] + 1 < lms_count)
                ? lms_positions[lms_map[current_lms] + 1]
                : length;

            bool same_substring = true;

            if (previous_end - previous_lms !=
                current_end - current_lms) {
                same_substring = false;
            } else {
                while (previous_lms < previous_end) {
                    if (text[previous_lms] != text[current_lms]) {
                        same_substring = false;
                        break;
                    }

                    previous_lms++;
                    current_lms++;
                }

                if (previous_lms == length ||
                    text[previous_lms] != text[current_lms]) {
                    same_substring = false;
                }
            }

            if (!same_substring)
                reduced_alphabet_max++;

            reduced_text[lms_map[sorted_lms_positions[i]]] =
                reduced_alphabet_max;
        }

        vector<int> reduced_suffix_array =
            sa_is(reduced_text, reduced_alphabet_max);

        for (int i = 0; i < lms_count; i++) {
            sorted_lms_positions[i] =
                lms_positions[reduced_suffix_array[i]];
        }

        induce(sorted_lms_positions);
    }

    return suffix_array;
}


int get_char(char nucleotide) {
    if (nucleotide == 'A') return 1;
    if (nucleotide == 'C') return 2;
    if (nucleotide == 'G') return 3;
    if (nucleotide == 'T') return 4;

    return 5;
}

char complement(char nucleotide) {
    if (nucleotide == 'A') return 'T';
    if (nucleotide == 'T') return 'A';
    if (nucleotide == 'C') return 'G';
    if (nucleotide == 'G') return 'C';

    return 'N';
}

string reverse_complement(const string& sequence) {
    string result;
    result.reserve(sequence.size());

    for (int i = (int)sequence.size() - 1; i >= 0; i--)
        result.push_back(complement(sequence[i]));

    return result;
}

string read_line(gzFile file) {
    string result;
    int character;

    while ((character = gzgetc(file)) != -1) {
        if (character == '\n')
            break;

        if (character != '\r')
            result.push_back((char)character);
    }

    return result;
}

string read_genome(const string& filename) {
    gzFile file = gzopen(filename.c_str(), "rb");

    if (!file)
        throw runtime_error("Cannot open genome file");

    string genome;
    string line;

    while (true) {
        line = read_line(file);

        if (line.empty()) {
            if (gzeof(file))
                break;

            continue;
        }

        if (line[0] == '>')
            continue;

        for (char nucleotide : line) {
            nucleotide = toupper((unsigned char)nucleotide);

            if (nucleotide == 'A' || nucleotide == 'C' ||
                nucleotide == 'G' || nucleotide == 'T' ||
                nucleotide == 'N') {
                genome.push_back(nucleotide);
            }
        }

        if (gzeof(file))
            break;
    }

    gzclose(file);

    return genome;
}

int compare_suffix(
    const string& genome,
    int suffix_position,
    const string& pattern
) {
    for (int i = 0; i < (int)pattern.size(); i++) {
        if (suffix_position + i >= (int)genome.size())
            return -1;

        int genome_character = get_char(genome[suffix_position + i]);
        int pattern_character = get_char(pattern[i]);

        if (genome_character < pattern_character)
            return -1;

        if (genome_character > pattern_character)
            return 1;
    }

    return 0;
}

void find_occurrences(
    const string& genome,
    const vector<int>& suffix_array,
    const string& pattern,
    vector<int>& positions
) {
    int suffix_array_size = suffix_array.size();

    int left_boundary = 0;
    int right_boundary = suffix_array_size;

    while (left_boundary < right_boundary) {
        int middle = left_boundary +
            (right_boundary - left_boundary) / 2;

        if (compare_suffix(
                genome, suffix_array[middle], pattern
            ) < 0) {
            left_boundary = middle + 1;
        } else {
            right_boundary = middle;
        }
    }

    int first_match = left_boundary;

    left_boundary = 0;
    right_boundary = suffix_array_size;

    while (left_boundary < right_boundary) {
        int middle = left_boundary +
            (right_boundary - left_boundary) / 2;

        if (compare_suffix(
                genome, suffix_array[middle], pattern
            ) <= 0) {
            left_boundary = middle + 1;
        } else {
            right_boundary = middle;
        }
    }

    int after_last_match = left_boundary;

    for (int i = first_match; i < after_last_match; i++) {
        int genome_position = suffix_array[i];

        if (genome_position + (int)pattern.size() <=
            (int)genome.size()) {
            positions.push_back(genome_position);
        }
    }
}

int main(int argc, char* argv[]) {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 3) {
        cerr << "Usage: ./mapper genome.fna.gz reads.fastq.gz\n";
        return 1;
    }

    string genome_filename = argv[1];
    string reads_filename = argv[2];

    cerr << "Reading reference genome...\n";

    string genome = read_genome(genome_filename);
    int genome_length = genome.size();

    cerr << "Genome length: " << genome_length << '\n';

    vector<int> encoded_genome;
    encoded_genome.reserve(genome_length + 1);

    for (char nucleotide : genome)
        encoded_genome.push_back(get_char(nucleotide));

    encoded_genome.push_back(0);

    cerr << "Building suffix array using SA-IS...\n";

    vector<int> suffix_array = sa_is(encoded_genome, 5);

    cerr << "Suffix array built.\n";

    encoded_genome.clear();
    encoded_genome.shrink_to_fit();

    vector<long long> coverage_difference_all(genome_length + 1, 0);
    vector<long long> coverage_difference_unique(genome_length + 1, 0);

    long long total_reads = 0;
    long long mapped_reads = 0;
    long long unique_reads = 0;
    long long multi_reads = 0;
    long long unmapped_reads = 0;

    long long total_placements = 0;
    long long total_aligned_bases = 0;

    cerr << "Processing FASTQ reads...\n";

    gzFile reads_file = gzopen(reads_filename.c_str(), "rb");

    if (!reads_file)
        throw runtime_error("Cannot open FASTQ file");

    while (true) {
        string header = read_line(reads_file);

        if (header.empty()) {
            if (gzeof(reads_file))
                break;

            continue;
        }

        string sequence = read_line(reads_file);
        string separator = read_line(reads_file);
        string quality = read_line(reads_file);

        if (header[0] != '@' ||
            separator.empty() || separator[0] != '+') {
            throw runtime_error("Invalid FASTQ record");
        }

        total_reads++;

        for (char& nucleotide : sequence)
            nucleotide = toupper((unsigned char)nucleotide);

        bool valid_sequence = true;

        for (char nucleotide : sequence) {
            if (nucleotide != 'A' && nucleotide != 'C' &&
                nucleotide != 'G' && nucleotide != 'T') {
                valid_sequence = false;
                break;
            }
        }

        if (!valid_sequence || sequence.empty() ||
            sequence.size() > genome.size()) {
            unmapped_reads++;
            continue;
        }

        vector<int> matched_positions;

        find_occurrences(
            genome,
            suffix_array,
            sequence,
            matched_positions
        );


        string reverse_complement_sequence =
            reverse_complement(sequence);

        if (reverse_complement_sequence != sequence) {
            find_occurrences(
                genome,
                suffix_array,
                reverse_complement_sequence,
                matched_positions
            );
        }

        sort(matched_positions.begin(), matched_positions.end());

        matched_positions.erase(
            unique(matched_positions.begin(), matched_positions.end()),
            matched_positions.end()
        );

        if (matched_positions.empty()) {
            unmapped_reads++;
            continue;
        }

        mapped_reads++;

        int read_length = sequence.size();

        if (matched_positions.size() == 1)
            unique_reads++;
        else
            multi_reads++;

        total_placements += matched_positions.size();

        // Update coverage for both uniquely mapped and multiply mapped.
        for (int genome_position : matched_positions) {
            coverage_difference_all[genome_position]++;
            coverage_difference_all[genome_position + read_length]--;

            total_aligned_bases += read_length;
        }

        // Update coverage only for uniquely mapped reads.
        if (matched_positions.size() == 1) {
            int genome_position = matched_positions[0];

            coverage_difference_unique[genome_position]++;
            coverage_difference_unique[
                genome_position + read_length
            ]--;
        }

        if (total_reads % 100000 == 0) {
            cerr << "Processed " << total_reads
                 << " reads. Mapped: " << mapped_reads
                 << '\n';
        }
    }

    gzclose(reads_file);

    long long covered_bases_all = 0;
    long long covered_bases_unique = 0;

    long long depth_sum_all = 0;
    long long depth_sum_unique = 0;

    long long current_coverage_all = 0;
    long long current_coverage_unique = 0;

    for (int i = 0; i < genome_length; i++) {
        current_coverage_all += coverage_difference_all[i];
        current_coverage_unique += coverage_difference_unique[i];

        depth_sum_all += current_coverage_all;
        depth_sum_unique += current_coverage_unique;

        if (current_coverage_all > 0)
            covered_bases_all++;

        if (current_coverage_unique > 0)
            covered_bases_unique++;
    }

    cout << fixed << setprecision(4);

    cout << "\nMAPPING REPORT\n";

    cout << "Genome length: "
         << genome_length << '\n';

    cout << "Total reads: "
         << total_reads << '\n';

    cout << "Mapped reads: "
         << mapped_reads << '\n';

    cout << "Unmapped reads: "
         << unmapped_reads << '\n';

    cout << "Unique reads: "
         << unique_reads << '\n';

    cout << "Multi-mapped reads: "
         << multi_reads << '\n';

    cout << "\nMapping rate: "
         << 100.0 * mapped_reads / max(1LL, total_reads)
         << "%\n";

    cout << "Unique mapping rate: "
         << 100.0 * unique_reads / max(1LL, total_reads)
         << "%\n";

    cout << "Multi-mapping rate: "
         << 100.0 * multi_reads / max(1LL, total_reads)
         << "%\n";

    cout << "Unmapped rate: "
         << 100.0 * unmapped_reads / max(1LL, total_reads)
         << "%\n";

    cout << "\nTotal alignment placements: "
         << total_placements << '\n';

    cout << "Genome coverage (all placements): "
         << 100.0 * covered_bases_all / genome_length << "%\n";

    cout << "Genome coverage (unique reads only): "
         << 100.0 * covered_bases_unique / genome_length << "%\n";

    cout << "\nMean aligned bases per mapped read: "
         << (double)total_aligned_bases /
            max(1LL, mapped_reads) << '\n';

    return 0;
}
