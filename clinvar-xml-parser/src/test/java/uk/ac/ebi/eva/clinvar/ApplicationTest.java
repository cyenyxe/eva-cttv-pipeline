/*
 *  Copyright 2018 EMBL - European Bioinformatics Institute
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
package uk.ac.ebi.eva.clinvar;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import static org.junit.Assert.assertEquals;

public class ApplicationTest {

    @Rule
    public TemporaryFolder outputFolder = new TemporaryFolder();

    @Test
    public void clinvarIntegrationTest() throws Exception {
        // application command line arguments
        Path inputFilePath = Paths.get(this.getClass().getResource("/ClinvarExample.xml.gz").toURI());
        String outputFolderAbsolutePath = outputFolder.getRoot().getAbsolutePath();
        String[] args = {"--inputFileName", inputFilePath.toString(),
                         "--outputDir", outputFolderAbsolutePath};

        // execute application
        Application.main(args);

        // check that an output gzipped file containing 10 lines has been generated
        GZIPInputStream gzipJsonFile = new GZIPInputStream(
                new FileInputStream(outputFolderAbsolutePath + "/" + Application.OUTPUT_FILE_NAME));
        BufferedReader br = new BufferedReader(new InputStreamReader(gzipJsonFile));
        assertEquals(10, br.lines().count());
    }

}